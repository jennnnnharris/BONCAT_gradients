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
  --i-data 16S.qza \
  --o-visualization 16S_raw.qzv



