#!/bin/bash
#SBATCH --account=one
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=6:00:00
#SBATCH --mem=36gb
#SBATCH --job-name=filter.tre
#SBATCH --output="filt.tr.out"
#SBATCH --error="filt.t.err"
#Get Started

echo "Job started on $(hostname) at $(date)"

#start conda
module load anaconda3
source activate qiime2-env

# go the right place
cd /storage/group/ltb5167/default/JennHarris/BONCAT_16S/qiime_sepp

qiime phylogeny filter-tree \
 --i-tree insertion-tree.qza \
 --i-table filtered_table.qza \
 --o-filtered-tree filter-tree.qza