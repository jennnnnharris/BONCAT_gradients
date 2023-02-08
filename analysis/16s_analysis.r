# 16s analyis

### 1. Initial Setup ###

### Set the working directory; modify to your own ###

setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")

### Clear workspace ###

rm(list=ls())

## if trouble loading phyloseq, see: http://joey711.github.io/phyloseq/install

#source("http://bioconductor.org/biocLite.R")
#biocLite("Heatplus")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

### Load required libraries ###

library(ade4)
library(vegan)
library(RColorBrewer)
library(phyloseq)
library(gplots)
library(Heatplus)
library(dplyr)
library(tidyverse)
library(viridis)
library(hrbrthemes)

### Import Data ###

taxon <- read.table("16s/taxonomy.tsv", sep="\t", header=T, row.names=1)
otus <- read.table("16s/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- readxl::read_excel("metadata.xlsx")


