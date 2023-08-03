########################differential abundance analysis #######################
# tutorialS
## from https://rpubs.com/lwaldron/EPIC2017_linearmodelinglab
## https://mcbl.readthedocs.io/en/master/tut-phyloseq.html
## https://www.yanh.org/2021/01/01/microbiome-r/#differential-abundance-analysis
## https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
## https://joey711.github.io/phyloseq-extensions/DESeq2.html
# jenn harris
# july 2023 
#install.packages("tidyverse")
library(tidyverse)
library(vegan)
library(phyloseq)
library(readxl)
library(DESeq2)
library("ggplot2")

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")


# Set the working directory; modify to your own ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
taxon <- read.table("asv_level_output/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)

##Data wrangleing for importing into Phyloseq##
#recode metadata 
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))
#Transpose asv table
asvs.t <- t(asvs.raw)
#order metadata
metadat<-metadat[order(metadat$SampleID),]
#order otu table
asvs.t<-asvs.t[order(row.names(asvs.t)),]
#order taxon table
taxon<-taxon[,1:7]
#make rownames sample names
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

##import non rarefied data into phyloseq ##
Workshop_OTU <- otu_table(asvs.t, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
#take a look at PS object
print(ps)
#we a get a Plyloseq object with  15027 taxa

##------ remove chloroplasts##
# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
# 14833 taxa
##----- remove taxa that have less than 10 count in all samples##
#This is a good way to improve your power to detect differentially abundant taxa 
#because these low-abundance taxa are less likely to be differentially abundant 
#and increase multiple testing.
ps<-prune_taxa(taxa_sums(ps) > 10, ps)
ps
sample_data(ps)
#8031 taxa and 42 samples 
head(sample_data(ps)$Fraction)

ps = subset_samples(ps, Compartment!= "ctl" & Fraction != "Total_DNA" & Compartment == "Rhizosphere")
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
head(sample_data(ps)$Fraction)
head(sample_data(ps)$Compartment)
ps
sample_data(ps)
##convert to DESeq2 object
dds.data = phyloseq_to_deseq2(ps, ~ Fraction)
dds.t = DESeq(dds.data, test="Wald", fitType="parametric")

res = results(dds.t, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

dim(sigtab)


windows(width = 6, height=7)

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$ Phyla, function(x) max(x))
x = sort(x, TRUE)
sigtab$ Phyla = factor(as.character(sigtab$ Phyla), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$ Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$ Family = factor(as.character(sigtab$ Family), levels=names(x))
ggplot(sigtab, aes(x= Family, y=log2FoldChange, color= Phyla)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

