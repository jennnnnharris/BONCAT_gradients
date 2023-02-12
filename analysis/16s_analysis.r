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
metadat <- read.delim("16s/metadata.txt", sep="\t", header = T, check.names=FALSE)


## Transpose OTU table ##
otus.t <- t(otus)

## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order otu table
otu.t<-otus.t[order(row.names(otus.t)),]

## Determine minimum available reads per sample ##

min(rowSums(otus.t))

### Rarefy to obtain even numbers of reads by sample ###

set.seed(336)
otus.r<-rrarefy(otus.t, 41610)

## Convert OTU numbers to percentages ##
otus.perc<-otus.r/rowSums(otus.r)*100


#------- 1. Run PCoA analysis of entire dataset ------

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otus.perc, method = "bray")

# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(nrow(otus.perc)-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #

otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

#to make coloring things easier I'm gong to added a combined fractionXboncat column no sure if i need this

metadat<-mutate(metadat, fraction_BCAT = paste0(metadat$Fraction, metadat$BONCAT))

## Plot ordination with factors coloured and shaped as you like # 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

svg(file="figures/16s/Pcoa.svg",width = 8, height=8 )
#windows(title="PCoA on OTUs - Bray Curtis")
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, col=c("black"),
       pch=c(21,21,22,23,24)[as.factor(metadat$Fraction)],
       lwd=1,cex=2,
       #bg=c("#003f5c","grey", "#bc5090", "#ffa600")[as.factor(metadat$BONCAT)])
       bg=c( "grey27","grey","#9795ff","#4406e2", "#ffb4f6","#e20a8f", "#695e00", "#6eda0a", "#0a7416" )[as.factor(metadat$fraction_BCAT)])

#dark brown, grey, endo (light purple, dark purple), nod(light pink, darkpink ), rhizo(darkest green, light green, darkgreen)
# BulkDNA ctlctl EndoPOS EndoSYBR NodPOS NodSYBR RhizoDNA RhizoPOS RhizoSYBR

legend("top",legend=c("BulkDNA", "ctl", "EndoPOS", "EndoSYBR", "NodPOS", "NodSYBR", "RhizoDNA", "RhizoPOS", "RhizoSYBR"),
  pch=c(16,16,15,15,18, 18,17,17,17),
  cex=1.1, 
  col=c("grey27","grey","#9795ff","#4406e2", "#ffb4f6","#e20a8f", "#695e00", "#6eda0a", "#0a7416" ))
dev.off()


# permanova

otu.perm<- adonis2(otus.perc~ Fraction*BONCAT, data = metadat, permutations = 999, method="bray")
otu.perm
# Fraction         4   8.3386 0.59802 18.4333  0.001 ***
#  BONCAT           2   1.2988 0.09315  5.7422  0.001 ***
#  Fraction:BONCAT  2   0.5742 0.04118  2.5387  0.126   

otu.ano<- anosim(otus.perc, grouping =  metadat$Fraction, permutations = 999)
summary(otu.ano)

dispersion <- betadisper(otus.bray, group=metadat$Fraction)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Groups     4 0.61336 0.153340 40.69    999  0.001 ***
 
# non homogeneous varience :( 

dispersion <- betadisper(otus.bray, group=metadat$BONCAT)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Groups     3 0.16064 0.053547 1.148    999  0.363

# homogenous varience :)

dispersion <- betadisper(otus.bray, group=metadat$fraction_BCAT)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Groups     8 0.69417 0.086771 26.697    999  0.001 ***
  
# non homogenous variance 



###### if I just look at rhizo is active verse no active different

test<-otus.perc[which(metadat$Fraction == "Rhizo"),]
metadat_t<-metadat[which(metadat$Fraction == "Rhizo"),]

otu.perm<- adonis2(test~ BONCAT, data = metadat_y, permutations = 999, method="bray")
otu.perm
#BONCAT    2  0.97631 0.40122 3.6854  0.001 ***

#nod
test<-otus.perc[which(metadat$Fraction == "Nod"),]
metadat_t<-metadat[which(metadat$Fraction == "Nod"),]

otu.perm<- adonis2(test~ BONCAT, data = metadat_t, permutations = 999, method="bray")
otu.perm
#BONCAT    1  0.11358 0.58017 9.6736  0.024 *

#nod
test<-otus.perc[which(metadat$Fraction == "Endo"),]
metadat_t<-metadat[which(metadat$Fraction == "Endo"),]

otu.perm<- adonis2(test~ BONCAT, data = metadat_t, permutations = 999, method="bray")
otu.perm
#BONCAT    1  0.19218 0.58686 9.9433  0.014 *
  

  

