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
qsub -I -A open -N filter -l nodes=1:ppn=10 -l pmem=8gb -l walltime=2:00:00
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

qiime cutadapt trim-paired \
    --i-demultiplexed-sequences 16S.qza  \
    --p-cores 1 \
    --p-adapter-f ^GTGYCAGCMGCCGCGGTAA...AAACTYAAAKRAATTGRCGG \
    --p-adapter-r ^CCGYCAATTYMTTTRAGTTT...TTACCGCGGCKGCTGRCAC \
    --p-error-rate 0.1 \
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

# copy to local device

#scp jeh6121@datamgr.aci.ics.psu.edu:/storage/work/jeh6121/TTC_16S/paired-end-demux.trimmed.primer.qzv .

# deniose with dada 2
# tutorial https://github.com/allenlab/QIIME2_16S_ASV_protocol

qiime dada2 denoise-paired \
    --i-demultiplexed-seqs paired-end-demux.trimmed.primer.qza \
    --p-n-threads 0 \
    --p-trunc-q 2 \
    --p-trunc-len-f 219 \
    --p-trunc-len-r 194 \
    --p-max-ee-f 2 \
    --p-max-ee-r 4 \
    --p-n-reads-learn 1000000 \
    --p-chimera-method pooled \
    --o-table table-dada2.qza \
    --o-representative-sequences demux.trimmed.dada2.qza \
    --o-denoising-stats stats-dada2.qza
# upload metadata
# download this file from teh cluster to view it
# scp jeh6121@datamgr.aci.ics.psu.edu:/storage/work/jeh6121/TTC_16S/paired-end-demux.trimmed.primer.qzv .

qiime tools export \
  --input-path demux.trimmed.dada2.qza \
  --output-path asv_table

qiime tools export \
  --input-path table-dada2.qza \
  --output-path asv_table

biom convert -i asv_table/feature-table.biom -o asv_table/feature-table.tsv --to-tsv

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization tabledada2.qzv \
  --m-sample-metadata-file metadata.txt

qiime feature-table tabulate-seqs \
  --i-data demux.trimmed.dada2.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

qiime tools export \
  --input-path demux.trimmed.dada2.qza \
  --output-path asvs

# silva classification 

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
  --i-classifier silva/silva-138-99-classifier.qza \
  --i-reads demux.trimmed.dada2.qza \
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


# phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences demux.trimmed.dada2.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


#download things from the cluster to my pc

# mkdir qzv_files
# mv *qsv qzv_files

scp -r jeh6121@datamgr.aci.ics.psu.edu:/storage/work/jeh6121/TTC_16S/asv_table .














