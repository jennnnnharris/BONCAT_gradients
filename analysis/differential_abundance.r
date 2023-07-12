########################differential abundance analysis #######################
# tutorial from https://rpubs.com/lwaldron/EPIC2017_linearmodelinglab
# jenn harris
# july 2023 
#install.packages("tidyverse")
library(tidyverse)
library(vegan)
library(phyloseq)
library(readxl)
install.packages("DESeq2")
library("DESeq2"); packageVersion("DESeq2")


# Set the working directory; modify to your own ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
taxon <- read.delim("otu_output/taxonomy.txt", sep="\t", header=T, row.names=1)
otu.raw <- read.table("otu_output/feature-table.txt", sep="\t", header = T , row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)

##Data wrangleing for importing into Phyloseq##
#recode metadata 
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))
#Transpose otu table
otu.t <- t(otu.raw)
#order metadata
metadat<-metadat[order(metadat$SampleID),]
#order otu table
otu.t<-otu.t[order(row.names(otu.t)),]
#order taxon table
taxon<-taxon[,1:7]
#make rownames sample names
metadat<-as.matrix(metadat)
y<-colnames(otu.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

##import non rarefied data into phyloseq ##
Workshop_OTU <- otu_table(otu.t, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
#take a look at PS object
print(ps)
#we a get a Plyloseq object with  7727 taxa

##------ remove chloroplasts##
# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)

##----- remove taxa that have less than 10 count in all samples##
#This is a good way to improve your power to detect differentially abundant taxa 
#because these low-abundance taxa are less likely to be differentially abundant 
#and increase multiple testing.
ps<-prune_taxa(taxa_sums(ps) > 10, ps)
ps
sample_data(ps)
#3934 taxa and 42 samples 

##convert to DESeq2 object
dds.data = phyloseq_to_deseq2(ps, ~country)
