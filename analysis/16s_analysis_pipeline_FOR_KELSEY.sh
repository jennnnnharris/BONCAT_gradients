#bioinformatics pipeline for 16s microbiome data

#you will need
## 16s sequences
## Qiime2

#step are
## 1. install software if you haven't already
## 2. make a qiime map file
## 3. check quality of reads
## 4. Pair reads and assign ASVs
## 5. Assign taxonomy. 




##### logon and install software ###########
#logon to PSU roar collab
ssh jeh6121@submit.hpc.psu.edu # change to your username
# load anaconda
module load anaconda3
# create conda environment
conda create -n qiime2-env -y
#start conda
module load anaconda3
# add wget
conda install wget
#load software
# These instructions are identical to the Linux (64-bit) instructions
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-env --file qiime2-2022.2-py38-linux-conda.yml
# get it's there
/storage/home/jeh6121/.conda/envs/qiime2-env

####### logging on and starting software ############
#logon to PSU roar collab
ssh jeh6121@submit.hpc.psu.edu # change to your username
#start interactive job
salloc -N 1 -n 4 --mem-per-cpu=1024 -t 2:00:00
# go to work directory
cd /storage/group/ltb5167/default/JennHarris/BONCAT_16S
#start conda
module load anaconda3
source activate qiime2-env
# if that doesn't work check what conda environemnt you already have. you should see qiime
cd /storage/home/jeh6121/.conda/envs

####### Make qiime mapping file ##########
# see example file for format
# upload to cluster

###### importing Sequences into Qiime
##### you will need a qiime mappping file for this.

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_mapping.txt \
  --output-path 16S.qza \
  --input-format PairedEndFastqManifestPhred33V2

#trimming adapters
# illumina AGATCGGAAGAGC

qiime cutadapt trim-paired \
  --p-cores 6 \
  --i-demultiplexed-sequences 16Strimmed.qza \
  --p-adapter-f GTGYCAGCMGCCGCGGTAA...GGACTACNVGGGTWTCTAA \
  --p-adapter-r GGACTACNVGGGTWTCTAA...TTAGAWCCCNGTAGTCC \
  --p-error-rate 0.1 \
  --o-trimmed-sequences 16S.trimmed.primer.qza

# view
qiime demux summarize \
    --i-data 16s.trimmed.primer.qza \
    --p-n 100000 \
    --o-visualization 16s.trimmed.primer.qzv
# view this file on qiimes website. trim where quality drops below 30 phred score.

# deniose with dada 2
# tutorial https://github.com/allenlab/QIIME2_16S_ASV_protocol

mkdir asvs

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ./16S.trimmed.primer.qza \
    --p-n-threads 0 \
    --p-trunc-q 2 \
    --p-trunc-len-f 205 \
    --p-trunc-len-r 205 \
    --p-max-ee-f 2 \
    --p-max-ee-r 5 \
    --p-n-reads-learn 1000000 \
    --p-chimera-method consensus \
    --o-table ./asvs/feature.table.dada2.qza \
    --o-representative-sequences ./asvs/rep.seqs.dada2.qza \
    --o-denoising-stats ./asvs/stats.dada2.qza

# Stats on how DADA2 went

qiime tools export \
   --input-path ./asvs/stats.dada2.qza \
   --output-path ./asvs/stats.dada2.qzv

# most samples had low chimeras, 80% to 90% of sequences were non chimeric
# for some samples (C1B DNA, C2B DNA, C2R DNA, C5B DNA) about 50% of sequences where chimeras
# this is not great but I'm going to keep moving forward the pipeline

#### summerize and export

qiime tools export \
  --input-path  table-dada2.qza\
  --output-path asvs/asv_table

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

# export representative sequences 
# 
qiime tools export \
  --input-path ./asvs/trimmed.dada2.qza\
  --output-path asvs

# tablulate and export feature table
qiime feature-table tabulate-seqs \
  --i-data  ./asvs/trimmed.dada2.qza \
  --o-visualization rep-seqs.qzv

#export with metadata
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization asvs/table-dada2.qzv \
  --m-sample-metadata-file metadata.tls
xt

qiime tools export \
  --input-path rep.seqs.dada2.qza\
  --output-path ./


######################## silva classification ############################

# downloaded the silva  classfier sequences from Qiimes website
# the upload files to cluster 
#scp silva-138-99-seqs-515-806.qza silva-138-99-tax-515-806.qza jeh6121@datamgr.aci.ics.psu.edu:/storage/work/jeh6121/TTC_16S/

# train classifier
# optimize clasifier by finding sections that are with the primers you used
qiime feature-classifier extract-reads \
  --i-sequences silva/silva-138-99-seqs-515-806.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer CCGYCAATTYMTTTRAGTTT \
  --o-reads silva/silva-138-99-extracts.qza

# train the classifiers
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /silva/silva-138-99-extracts.qza \
  --i-reference-taxonomy /silva/silva-138-99-tax.qza \
  --o-classifier /silva/silva-138-99-classifier.qza

# test the classifier
qiime feature-classifier classify-sklearn
  --i-classifier silva/silva-138-99-nb-classifier.qza \
  --i-reads otu/rep.seqs.otu.qza \
  --o-classification otu/taxonomy.otu2.qza

rep.seqs.otu.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv


######## using the pre trained classifier #######

qiime feature-classifier classify-sklearn \
 --i-classifier 2022.10.backbone.v4.nb.qza \
  --i-reads demux.trimmed.dada2.qza \
  --o-classification taxonomy.qza


#######################download things from the cluster to my pc#########################

mkdir qzv_files
mv *.qzv qzv_files


#from the cluster I need:
#1 metadata (usually called metadate.txt)
#2 table of all the asvs (usually feature-table.tsv)
#3 taxomony file (usually taxnomy.tsv)
#4 export representative sequences (called asvs/trimmed.dada2.tsv)

qiime tools export \
--input-path feature.table.otu.qza \
--output-path  /storage/group/ltb5167/default/JennHarris/BONCAT_16S/otu/output

qiime tools export \
--input-path table-dada2.qza \
--output-path /storage/group/ltb5167/default/JennHarris/BONCAT_16S/output_greengenes2




