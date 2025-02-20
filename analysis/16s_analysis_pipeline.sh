#bioinformatics pipeline for 16s microbiome data
#this script is in bash 
#designed to be run line by line in the terminal on linux

#you will need
## 16s sequences
## Qiime2

#step are
## 1. install software 
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
    --p-trunc-len-f 220 \
    --p-trunc-len-r 220 \
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
# for some samples (C1B DNA, C2B DNA, C2R DNA, C5B DNA) about 60% of total sequences where chimeras
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



############### green genes taxanomy classification #############################

pip install q2-greengenes2
wget http://ftp.microbio.me/greengenes_release/2022.10/2022.10.taxonomy.asv.nwk.qza

qiime greengenes2 filter-features \
	--i-feature-table feature-table.biom \
	--i-reference 2022.10.taxonomy.asv.nwk.qza \
	--o-filtered-feature-table filter.table.biom.qza

qiime greengenes2 taxonomy-from-table \
     --i-reference-taxonomy 2022.10.taxonomy.asv.nwk.qza \
     --i-table filter.table.biom.qza \
     --o-classification table.taxonomy.qza

wget http://ftp.microbio.me/greengenes_release/current/2022.10.phylogeny.asv.nwk.qza

# green genes2 has a tree I can probably use


######## using the pre trained classifier #######

qiime feature-classifier classify-sklearn \
 --i-classifier 2022.10.backbone.v4.nb.qza \
  --i-reads demux.trimmed.dada2.qza \
  --o-classification taxonomy.qza

####### assign phylogeny ##########

##use MAFFT to make a reference alignment

qiime alignment mafft \
--i-sequences 2022.10.backbone.v4.fna.qza \
--o-alignment aligned-ref-sequences.qza

# another helpful tutorial https://john-quensen.com/tutorials/processing-16s-sequences-with-qiime2-and-dada2/
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences asvs/trimmed.dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime tools export \
--input-path unrooted.tree.qza \
--output-path /storage/group/ltb5167/default/JennHarris/BONCAT_16S/otu/output


qiime tools export \
--input-path rooted-tree.qza \
--output-path /storage/group/ltb5167/default/JennHarris/BONCAT_16S/otu/output

#######################download things from the cluster to my pc#########################

mkdir qzv_files
mv *.qzv qzv_files


#from the cluster I need:
#1 metadata (usually called metadate.txt)
#2 table of all the asvs (usually feature-table.tsv)
#3 taxomony file (usually taxnomy.tsv)

#tree
#4 phylogeny tree output rooted tree
#5 phylogeny unrooted tree
#6 export representative sequences (called asvs/trimmed.dada2.tsv)

qiime tools export \
--input-path feature.table.otu.qza \
--output-path  /storage/group/ltb5167/default/JennHarris/BONCAT_16S/otu/output

qiime tools export \
--input-path table-dada2.qza \
--output-path /storage/group/ltb5167/default/JennHarris/BONCAT_16S/output_greengenes2




