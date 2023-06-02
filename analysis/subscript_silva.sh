#!/bin/bash

#SBATCH -A open
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=12:00:00
#SBATCH --mem=12gb
#SBATCH --job-name=SILVA

#Get Started
echo "Job started on $(hostname) at $(date)"

module purge
module load anaconda3
source activate qiime2-env

cd /storage/group/ltb5167/default/JennHarris/criticalzone16Ssequencing_Jan2023/

# optimize clasifier by finding sections that are with the primers you used
#qiime feature-classifier extract-reads \
#  --i-sequences silva/silva-138-99-seqs-515-806.qza \
#  --p-f-primer GTGYCAGCMGCCGCGGTAA \
#  --p-r-primer CCGYCAATTYMTTTRAGTTT \
#  --o-reads silva/silva-138-99-extracts.qza

# train the classifiers
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads /silva/silva-138-99-extracts.qza \
#  --i-reference-taxonomy /silva/silva-138-99-tax.qza \
#  --o-classifier /silva/silva-138-99-classifier.qza


# test the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier /storage/group/ltb5167/default/JennHarris/silva/silva-138-99-nb-classifier.qza \
  --i-reads ./asvs/rep.seqs.dada2.qza \
  --o-classification ./asvs/taxonomy.qza

qiime metadata tabulate \
  --m-input-file ./asvs/taxonomy.qza \
  --o-visualization ./asvs/taxonomy.qzv


echo "finished"