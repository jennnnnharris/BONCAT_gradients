#bioinformatics pipeline for 16s microbiome data

#you will need
## 16s sequences
## Qiime2

#step are
## 1. install software if you haven't already
## 2. make a qiime map file
## 3. check the number of reads for each sample

####### logging on and starting software ############
#logon to PSU cluster
ssh jeh6121@submit.aci.ics.psu.edu
#start interactive job
qsub -I -A open -N filter -l nodes=1:ppn=10 -l pmem=8gb -l walltime=1:00:00
# go to work directory
cd /gpfs/group/ltb5167/default/JennHarris/BONCAT_16S
#start conda
module load anaconda3
source activate qiime2-2022.2
# if that doesn't work check what conda environemnt you already have. I couldn't remember if I install qiime
cd /storage/home/jeh6121/.conda/envs

####### Make qiime mapping file ##########
# see example file for format
# upload to cluster

###### importing Sequences into Qiime

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path qiime_mapping.txt \
  --output-path 16S.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data 16S.trimmed.primer.qza \
  --o-visualization 16S_trimmed.qzv


#trimming adapters
# illumina AGATCGGAAGAGC

# slightly edited:

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences 16S.qza  \
    --p-cores -6 \
    --p-adapter-f GTGYCAGCMGCCGCGGTAA \
    --p-adapter-r GGACTACNVGGGTWTCTAA \
    --p-error-rate 0.2 \
    --p-match-adapter-wildcards \
    --p-discard-untrimmed \
    --p-overlap 3 \
    --verbose \
    --o-trimmed-sequences 16Strimmed.qza

# tutorial here https://xyz1396.github.io/2021/09/26/16S-Tutorial/

# remove primers

#515f Modified	GTGYCAGCMGCCGCGGTAA	Parada et 
#806r Modified	GGACTACNVGGGTWTCTAA	Aprille
# reseve compliment of 806 GGACTACNVGGGTWTCTAA
#TTACCGCGGCGCTGCAC
#library(Biostrings)
# CCGTCAATTCMTTTRAGTTT...CCGCGGCKGCTGGCAC
#reverseComplement(DNAString("GTGCCAGCMGCCGCGG...AAACTYAAAKGAATTGACGG"))

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
    --p-max-ee-r 4 \
    --p-n-reads-learn 1000000 \
    --p-chimera-method consensus \
    --o-table ./asvs/table-dada2.qza \
    --o-representative-sequences ./asvs/trimmed.dada2.qza \
    --o-denoising-stats ./asvs/stats.dada2.qza

# on how DADA2 went

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

biom convert -i asvs/asv_table/feature-table.biom -o asvs/asv_table/feature-table.tsv --to-tsv

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
  --m-sample-metadata-file metadata.txt


######################## silva classification ############################

#upload files to cluster 
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
qiime feature-classifier classify-sklearn \
  --i-classifier silva/silva-138-99-nb-classifier.qza \
  --i-reads asvs/trimmed.dada2.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# if classifying with silva didn't work i'll just use the naive bayes trained on greenhenes sequences
# https://docs.qiime2.org/2022.2/tutorials/moving-pictures/index.html
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2022.2/common/gg-13-8-99-515-806-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads demux.trimmed.dada2.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

#output taxomny file

qiime tools export \
  --input-path taxonomy.qza \
  --output-path ./


####### assign phylogeny ##########

# another helpful tutorial https://john-quensen.com/tutorials/processing-16s-sequences-with-qiime2-and-dada2/


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences asvs/trimmed.dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime tools export \
--input-path unrooted-tree.qza \
--output-path phyloseq
qiime tools export \
--input-path rooted-tree.qza \
--output-path phyloseq

#download things from the cluster to my pc

mkdir qzv_files
mv *.qzv qzv_files

#from the cluster I need:
#1 metadata (usually called metadate.txt)
#2 table of all the asvs (usually feature-table.tsv)
#3 taxomony file (usually taxnomy.tsv)

#optional
#4 phylogeny tree output rooted tree
#5 phylogeny unrooted tree
#6 export representative sequences (called asvs/trimmed.dada2.tsv)
qiime tools export \
--input-path trimmed.dada2.qza \
--output-path  /gpfs/group/ltb5167/default/JennHarris/BONCAT_16S/phyloseq












