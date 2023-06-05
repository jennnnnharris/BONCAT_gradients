# 16s analysis

### 1. Initial Setup ###

### Clear workspace ###

rm(list=ls())

## if trouble loading phyloseq, see: http://joey711.github.io/phyloseq/install

#source("http://bioconductor.org/biocLite.R")
#biocLite("Heatplus")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("Heatplus")


### Load required libraries ###

#install.packages("tidyverse")
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(scales)
library(data.table)
library(phyloseq)
library(DT)
library(Heatplus)
library(viridis)
library(hrbrthemes)
library(ade4)
#library (gplots)
library(ape)

############-----load functions-------------######

# load in the function for making a heatmap with the tree #
heatmap.phylo <- function(x, Rowp, Colp, ...) {
  l = length(seq(-9.5, 10, .5))
  pal = colorRampPalette(c("#0c3fb6" ,"#cfc3c8", "#bb0000"))(l)
  row_order = Rowp$tip.label[Rowp$edge[Rowp$edge[, 2] <= Ntip(Rowp), 2]] 
  col_order = Colp$tip.label[Colp$edge[Colp$edge[, 2] <= Ntip(Colp), 2]] 
  x <- x[row_order, col_order]
  xl <- c(0.5, ncol(x)+0.5)
  yl <- c(0.5, nrow(x)+0.5)
  layout(matrix(c(0,1,0, 2,3,4, 0,5,0), nrow=3, byrow=TRUE),
         width=c(3.5,    4.5, 1),
         height=c(0.2, 3, 0.18))
  par(mar=rep(0,4))
  plot(Colp, direction="downwards", show.tip.label=FALSE,
       xaxs="i", x.lim=xl)
  par(mar=rep(0,4))
  plot(Rowp, direction="rightwards", show.tip.label=FALSE, 
       yaxs="i", y.lim=yl)
  lpp = .PlotPhyloEnv$last_plot.phylo 
  segments(lpp$xx[1:Ntip(Rowp)], lpp$yy[1:Ntip(Rowp)], par('usr')[2],
           lpp$yy[1:Ntip(Rowp)], lty=3, col='grey50')
  par(mar=rep(0,4), xpd=TRUE)
  image((1:ncol(x))-0.5, (1:nrow(x))-0.5, t(x), col=pal,
        xaxs="i", yaxs="i", axes=FALSE, xlab= "", ylab= "", breaks=seq(-10,10,.5))
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", yaxs="i", xlim=c(0,2), ylim=yl)
  text(rep(0,nrow(x)),1:nrow(x), row_order, pos=4,
       cex=1, xpd=NA)
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", xaxs="i", ylim=c(0,2), xlim=xl)
  text(1:ncol(x),rep(2,ncol(x)), col_order, srt=90, adj=c(1,.5),
       cex=1.5)
}

######-------------import Asvs data -----------------##############
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
taxon <- read.table("asv_level_output/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)

## Transpose ASVS table ##
asvs.t <- t(asvs.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs.t<-asvs.t[order(row.names(asvs.t)),]

# rearrange asv table
asvs.phyloseq<- (asvs.t)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_ASVS <- otu_table(asvs.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps.a <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

# take a look at PS object
print(ps.a)
# we a get a Plyloseq object with  15027 taxa

## Determine minimum available reads per sample ##
min.s<-min(rowSums(asvs.t))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs.t, min.s)
############------------ rarefaction curve------------############
S <- specnumber(otus.t) # observed number of species
raremax <- min(rowSums(otus.t))
plot(otus.t, otus.r, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
out<-rarecurve(otus.t, step = 20, sample = raremax, col = "blue", cex = 0.6)
## build plots
col <- c("black")
lty <- c("solid")
lwd <- rep(1, 42)
pars <- expand.grid(col = col, lty = lty, lwd = lwd, 
                    stringsAsFactors = FALSE)
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)
# save
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/rarecurve.svg",width = 6, height=6 )
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Sequences",
     ylab = "ASV", type = "n")
abline(v = raremax)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}
dev.off()


######--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

#####------make phyloseq object with rarefied data -------#####

asvs.phyloseq<- (asvs.r)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_ASVS <- otu_table(asvs.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps.a <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

#test it worked
#sample_names(ps)
print(ps.a)
# 15027 taxa

# remove chloroplast DNA
ps.a<-subset_taxa(ps.a, Class!=" Chloroplast")
ps.a<-subset_taxa(ps.a, Genus!=" Mitochondria")
ps.a<-subset_taxa(ps.a, Genus!=" Chloroplast")
# get rid of taxa that aren; in any samples
ps.a<-prune_taxa(taxa_sums(ps.a) > 0, ps.a)
any(taxa_sums(ps.a) == 0)
ps.a
# 14593 taxa

asvs.clean<-as.data.frame(t(as.data.frame(otu_table(ps.a))))
taxon<-as.data.frame(tax_table(ps.a))



####### import OTU 97% data##################
  # I chose to cluster into OTUs because when I looked at the sequencing inside the nodules
  # there were many asvs that all mapped to the rhizobia genus.
  # Those are all likely 1 species, given I they all BLAST to rhizobium leguminosarium v. trifoilia with WGS
  # I want to reduce the organism in the alpha diversity that are due to the dada2 clustering methods
  
## Set the working directory; modify to your own ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
taxon <- read.delim("otu_output/taxonomy.txt", sep="\t", header=T, row.names=1)
otu.raw <- read.table("otu_output/feature-table.txt", sep="\t", header = T , row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)


#data wrangleing and importing into Phyloseq#
## recode metadata 
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

## Transpose ASVS table ##
otu.t <- t(otu.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
otu.t<-otu.t[order(row.names(otu.t)),]


########-----------import non rarefied data into phyloseq #
## rearrange data for phyloseq
otu.phyloseq<- (otu.t)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otu.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

## import it phyloseq
Workshop_OTU <- otu_table(otu.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
##  take a look at PS object
print(ps)
# we a get a Plyloseq object with  7727 taxa

#######------ remove chloroplasts#
# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
# 7623 taxa
# output df 
otus<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))

####################rarefying taxa and importing into phyloseq #####
## Determine minimum available reads per sample ##
min.seqs<-min(rowSums(otu.t))
## Rarefy to obtain even numbers of reads by sample ##
set.seed(336)
otu.r<-rrarefy(otu.t, 41351)
## rearrange data for phyloseq
otu.phyloseq<- (otu.r)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otu.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)
## import it phyloseq
Workshop_OTU <- otu_table(otu.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps.r <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
ps.r
# we a get a Plyloseq object with  7727 taxa
#######------ remove chloroplasts-------#
# remove chloroplast DNA
ps.r<-subset_taxa(ps.r, Class!=" Chloroplast")
ps.r<-subset_taxa(ps.r, Genus!=" Mitochondria")
ps.r<-subset_taxa(ps.r, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
ps.r<-prune_taxa(taxa_sums(ps.r) > 0, ps.r)
any(taxa_sums(ps.r) == 0)
ps.r

# 7516 taxa
# output df 
otu.r.clean<-as.data.frame(t(as.data.frame(otu_table(ps.r))))
#taxon<-as.data.frame(tax_table(ps))

#####################------normalizing data by number of reads#################

## Set the working directory; modify to your own ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
taxon <- read.delim("otu_output/taxonomy.txt", sep="\t", header=T, row.names=1)
otu.raw <- read.table("otu_output/feature-table.txt", sep="\t", header = T , row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)


#data wrangleing and importing into Phyloseq#
## recode metadata 
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

## Transpose ASVS table ##
otu.t <- t(otu.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
otu.t<-otu.t[order(row.names(otu.t)),]


# hellinger data transformation 
# divides each value by the row sum and takes the square root of the quotient
# this would normalize by the number of reads in each sample and make very not abundant taxa less strongly influence the data
otu.n<-decostand(otu.t, "hellinger")

## rearrange data for phyloseq
otu.phyloseq<- (otu.n)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otu.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)
## import it phyloseq
Workshop_OTU <- otu_table(otu.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
ps
# 7727 taxa
#######------ remove chloroplasts-------#
# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps

# 7623 taxa
# output df 
otu.clean<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))






########------diversity figure---------########
# set wd for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

# diversity 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

svg(file="figures/16s/diversity.svg",width = 10, height=4 )
windows()
plot_richness(ps, "Fraction", measures = c("Observed","Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(fill = "Fraction")) + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578"))
dev.off()

# Data wrangling fo rdiversity of active microbes in each fraction
rich<-cbind(rich, metadat)
rich<-as.data.frame(rich)
colnames(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))
reorder(compartment_BCAT, -Observed)

rich<- rich %>%  filter(Plant!="NOPLANT", Fraction!="Inactive")
rich$compartment_BCAT <-factor(rich$compartment_BCAT, levels = c("Bulk_SoilTotal_DNA", "RhizosphereTotal_DNA", "RhizosphereTotal_Cells", "RhizosphereBONCAT_Active",
                                                                 "RootsTotal_Cells"   ,  "RootsBONCAT_Active" , "NoduleTotal_Cells" ,  "NoduleBONCAT_Active" ))          

# total + active number otus
svg(file="figures/16s/observed_divserity.svg",width = 6, height=4 )
windows(width = 10, height=7)
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=Observed, fill = Fraction, col= Fraction))+
  geom_boxplot() +
  scale_colour_manual(values = c( "orange",  "black", "black"))+
  scale_fill_manual( values = c("gold", "grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14), axis.text.y = element_text(size = 14) )+
  ylab("N Otus")+
  xlab("Compartment")
dev.off()



# total + active shannon
svg(file="figures/16s/shannon_diversity.svg",width = 6, height=4 )
windows(width = 10, height=7)
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=Shannon, fill = Fraction, col= Fraction))+
  geom_boxplot() +
  scale_colour_manual(values = c( "orange",  "black", "black"))+
  scale_fill_manual( values = c("gold", "grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14), axis.text.y = element_text(size = 14) )+
  ylab("shannon diversity")+
  xlab("Compartment")
dev.off()



# total + active simpson
svg(file="figures/16s/simpson_diversity.svg",width = 6, height=4 )
windows(width = 10, height=7)
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=Simpson, fill = Fraction, col= Fraction))+
  geom_boxplot() +
  scale_colour_manual(values = c( "orange",  "black", "black"))+
  scale_fill_manual( values = c("gold", "grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14), axis.text.y = element_text(size = 14) )+
  ylab("Simpson diversity")+
  xlab("Compartment")
dev.off()



############ taxa exploration in  the endophyte ######################

rich %>%
  filter(BONCAT!="DNA", Fraction!="ctl", Compartment== "Roots", REP!=1)%>%
  ggplot(aes(x=Fraction, y=Observed, col=REP))+
  #geom_boxplot() +
  #geom_line( x=Fraction, y=Observed, col= rep )
  #scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  ylab("Number of OTUS")+
  theme_bw()


rich %>%
  filter(BONCAT!="DNA", Fraction!="ctl", Compartment== "Roots")%>%
  ggplot(aes(x=Fraction, y=Shannon, col=REP))+
  #geom_boxplot() +
  #geom_line( x=Fraction, y=Observed, col= rep )
  #scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  ylab("Shannon divsersity")+
  theme_bw()


rich %>%
  filter(BONCAT!="DNA", Fraction!="ctl", Compartment== "Nodule")%>%
  ggplot(aes(x=Fraction, y=Shannon, col=REP))+
  #geom_boxplot() +
  #geom_line( x=Fraction, y=Observed, col= rep )
  #scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  ylab("Shannon divsersity")+
  theme_bw()


root_total<-subset_samples(ps, compartment_BCAT=="RootsTotal_Cells")
root_total<-prune_taxa(taxa_sums(root_total) > 0, root_total)
any(taxa_sums(root_total) == 0)
root_total
ot<-as.data.table(otu_table(root_total))
tt<-as.data.table(tax_table(root_total))
ot<-t(ot)
total<-cbind(ot,tt)
sum<-rowSums(total[,1:4])
total<-mutate(total, sum= sum)

window("total")
hist(total$sum)

root_active<-subset_samples(ps, compartment_BCAT=="RootsBONCAT_Active")
root_active<-prune_taxa(taxa_sums(root_active) > 0, root_active)
any(taxa_sums(root_active) == 0)
root_active
tt<-as.data.table(tax_table(root_active))
ot<-as.data.table(otu_table(root_active))
ot<-t(ot)
root_active<-cbind(ot,tt)
sum<-rowSums(root_active[,1:5])

root_active<-mutate(root_active, sum= sum)

windows()
hist(root_active$sum)
dim(total)
dim(root_active)
dim(filter(total, sum>1))
dim(filter(root_active, sum>1))

#What;s in the nodule?
nod<-subset_samples(ps, Compartment=="Nodule")
nod<-prune_taxa(taxa_sums(nod) > 0, nod)
any(taxa_sums(nod) == 0)
nod
# 40 taxa
ot<-as.data.table(otu_table(nod_active))
tt<-as.data.frame(tax_table(nod_active))
ot<-t(ot)
nod_active<-cbind(ot,tt)
sum<-rowSums(nod_active[,1:4])
nod_active<-mutate(nod_active, sum= sum)
# most everything in rhizosbia
window("total")
hist(nod_active$sum)


############################PCOA plots ########################
#Pcoa on rarefied Otu Data

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps.r), method = "bray")
ps.r


# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(42-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

#plant verse soil
metadat
# subset

#p<-otus.p[which(metadat$BONCAT != "ctl",)]

# subset metadata
as.factor(metadat$Compartment)
as.factor(metadat$Fraction)
as.factor(metadat$compartment_BCAT)
levels(as.factor(metadat$compartment_BCAT))
levels(as.factor(metadat$Fraction))

# colors :)
gold <- "#FFB000"
purple <- "#785EF0"
blue <- "#739AFF"
pink <- "#DC267F"

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/Pcoa_plantvsoil_raw.svg",width = 7, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, col=c("black"),
       pch=c(22, 21,22,24, 23)[as.factor(metadat$Fraction)],
       lwd=1,cex=2,
       bg=c( blue, "black", "white", pink , pink, purple, purple, purple, gold, gold)[as.factor(metadat$compartment_BCAT)])


legend("top",legend=c("beads", "Bulk_SoilTotal_DNA" ,      "neg control"   ,  "NoduleBONCAT_Active"    ,  "Nodule Total Cells"     ,      "RhizosphereBONCAT_Active",
                      "Rhizosphere Total Cells","RhizosphereTotal_DNA" ,    "Roots BONCAT Active"   ,    "Roots Total Cells"),
       pch=c(15,5,0,1,2 , 1, 2 , 5, 1, 2),
       col= c("black", "#739AFF", "black", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000"),
       bty = "n",
       inset = c(.05, 0))

dev.off()



##########-------run a new pcoa on just soil <3
ps.r
ps2<-subset_samples(ps.r, Compartment !=  "Nodule" & Compartment != "Roots")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample.names(ps2)
# 7401 taXA
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")

# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(24-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots") 

as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

square <- 22
diamond <- 23
triangle <- 24
circle <- 21

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/soil_raw.svg",width = 6, height=6 )
windows(title="PCoA on asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(square, circle, square, triangle, diamond)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c(blue,"black", "white", purple , purple, purple)[as.factor(metadat2$compartment_BCAT)])

legend("topleft",legend=c("Flow cyto control", "Bulk soil total DNA", " PCR control",  "Rhizosphere BONCAT_Active" , "Rhizosphere Inactive",  "Rhizosphere Total DNA"), 
       pch=c(15,5, 0, 1,2,5),
       cex=1.1, 
       col=c("black", "#739AFF",  "black", "#785EF0", "#785EF0",  "#785EF0"),
       bty = "n")

dev.off()

################# just active verse total

ps.r
ps2<-subset_samples(ps.r, Compartment !=  "Nodule" & Compartment != "Roots" & Fraction!="Total_DNA")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample.names(ps2)
# 7401 taXA
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")

# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(11-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots" & Fraction!="Total_DNA") 

as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

square <- 22
diamond <- 23
triangle <- 24
circle <- 21

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/soil_simple_raw.svg",width = 6, height=6 )
windows(title="PCoA on asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(square, circle, square, triangle)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c("black", "white", purple, purple)[as.factor(metadat2$compartment_BCAT)])

legend("topleft",legend=c("Flow cyto control", "Bulk soil total DNA", " PCR control",  "Rhizosphere BONCAT_Active" , "Rhizosphere Inactive",  "Rhizosphere Total DNA"), 
       pch=c(15,5, 0, 1,2,5),
       cex=1.1, 
       col=c("black", "#739AFF",  "black", "#785EF0", "#785EF0",  "#785EF0"),
       bty = "n")

dev.off()
#################-------------------just plant##############
ps.r
ps2<-subset_samples(ps.r, Compartment !=  "Bulk_Soil" & Compartment != "Rhizosphere" & Compartment != "ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample.names(ps2)
# 446 taXA
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")

# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(15-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Bulk_Soil" & Compartment != "Rhizosphere" & Compartment != "ctl") 

as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

square <- 22
diamond <- 23
triangle <- 24
circle <- 21

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/plant_raw.svg",width = 6, height=6 )
windows(title="PCoA on asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(circle, triangle)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c(pink, pink, gold, gold)[as.factor(metadat2$compartment_BCAT)])

legend("topleft",legend=c("Flow cyto control", "Bulk soil total DNA", " PCR control",  "Rhizosphere BONCAT_Active" , "Rhizosphere Inactive",  "Rhizosphere Total DNA"), 
       pch=c(15,5, 0, 1,2,5),
       cex=1.1, 
       col=c("black", "#739AFF",  "black", "#785EF0", "#785EF0",  "#785EF0"),
       bty = "n")

dev.off()


#########-------------- permanova----------------#########

# use output that doesn't have chloroplast for permanova. 


otu.r.clean.t <- t(otu.r.clean)
otu.perm<- adonis2(t(otu.r.clean)~ Compartment*Fraction, data = metadat, permutations = 999, method="bray")

otu.perm
# Fraction         4   8.3386 0.59802 18.4333  0.001 ***
#  BONCAT           2   1.2988 0.09315  5.7422  0.001 ***
#  Fraction:BONCAT  2   0.5742 0.04118  2.5387  0.126   

#analysis of similarities
otu.ano<- anosim(t(otu.r.clean), grouping =  metadat$Compartment, permutations = 999)
summary(otu.ano)

#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=metadat$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse



##### subset data by fraction 
#rhizo active verse total cells
test<-otu.r.clean.t[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_DNA"),]
metadat_y<-metadat[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_DNA") ,]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm

#rhizo active ver total DNA
test<-otu.r.clean.t[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive"),]
metadat_y<-metadat[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive") ,]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm


#nod active verse inactive
test<-otu.p.t[which(metadat$Compartment == "Nodule" &  metadat$Fraction != "Total_Cells" ),]
metadat_y<-metadat[which(metadat$Compartment == "Nodule" & metadat$Fraction != "Total_Cells" ),]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm

#endo active verse inactive
test<-otu.p.t[which(metadat$Compartment== "Roots" & metadat$Fraction != "Total_Cells" ),]
metadat_y<-metadat[which(metadat$Compartment == "Roots" & metadat$Fraction != "Total_Cells" ),]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm


#active roots verse nodules
test<-otu.p.t[which(metadat$Compartment != "Rhizosphere" &  metadat$Compartment != "Bulk_Soil" & metadat$Fraction == "BONCAT_Active"),]
metadat_y<-metadat[which(metadat$Compartment != "Rhizosphere" & metadat$Compartment != "Bulk_Soil" &   metadat$Fraction == "BONCAT_Active") ,]
otu.perm<- adonis2(test~ Compartment, data = metadat_y, permutations = 999, method="bray")
otu.perm

#active roots verse rhizo
test<-otu.p.t[which(metadat$Compartment != "Rhizosphere" &  metadat$Compartment != "Bulk_Soil" & metadat$Fraction == "BONCAT_Active"),]
metadat_y<-metadat[which(metadat$Compartment != "Rhizosphere" & metadat$Compartment != "Bulk_Soil" &   metadat$Fraction == "BONCAT_Active") ,]
otu.perm<- adonis2(test~ Compartment, data = metadat_y, permutations = 999, method="bray")
otu.perm


##################### PCOA ASVS level data  ##############

#Pcoa on rarefied ASVS Data

# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(otu_table(ps.a), method = "bray")
ps.a


# Perform PCoA analysis of BC distances #
asvs.pcoa <- cmdscale(asvs.bray, k=(42-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
asvs.p <- asvs.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
asvs.eig<-asvs.pcoa$eig
perc.exp<-asvs.eig/(sum(asvs.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

#plant verse soil
metadat
# subset

#p<-otus.p[which(metadat$BONCAT != "ctl",)]

# subset metadata
as.factor(metadat$Compartment)
as.factor(metadat$Fraction)
as.factor(metadat$compartment_BCAT)
levels(as.factor(metadat$compartment_BCAT))
levels(as.factor(metadat$Fraction))

# colors :)
gold <- "#FFB000"
purple <- "#785EF0"
blue <- "#739AFF"
pink <- "#DC267F"

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/Pcoa_ASVS_plantvsoil_raw.svg",width = 7, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(asvs.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(asvs.p, col=c("black"),
       pch=c(22, 21,22,24, 23)[as.factor(metadat$Fraction)],
       lwd=1,cex=2,
       bg=c( blue, "black", "white", pink , pink, purple, purple, purple, gold, gold)[as.factor(metadat$compartment_BCAT)])


legend("top",legend=c("beads", "Bulk_SoilTotal_DNA" ,      "neg control"   ,  "NoduleBONCAT_Active"    ,  "Nodule Total Cells"     ,      "RhizosphereBONCAT_Active",
                      "Rhizosphere Total Cells","RhizosphereTotal_DNA" ,    "Roots BONCAT Active"   ,    "Roots Total Cells"),
       pch=c(15,5,0,1,2 , 1, 2 , 5, 1, 2),
       col= c("black", "#739AFF", "black", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000"),
       bty = "n",
       inset = c(.05, 0))

dev.off()


##########-------run a new pcoa on just soil <3  ASVS level ###
ps.a
ps2a<-subset_samples(ps.a, Compartment !=  "Nodule" & Compartment != "Roots")
ps2a<-prune_taxa(taxa_sums(ps2a) > 0, ps2a)
any(taxa_sums(ps2a) == 0)
ps2a
sample.names(ps2a)
# 7401 taXA
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(otu_table(ps2a), method = "bray")

# Perform PCoA analysis of BC distances #
asvs.pcoa <- cmdscale(asvs.bray, k=(24-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
asvs.p <- asvs.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
asvs.eig<-asvs.pcoa$eig
perc.exp<-asvs.eig/(sum(asvs.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots") 

as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

square <- 22
diamond <- 23
triangle <- 24
circle <- 21

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/soil_raw.svg",width = 6, height=6 )
windows(title="PCoA on asvs- Bray Curtis", width = 7, height = 6)
ordiplot(asvs.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(asvs.p[,1:2],
       pch=c(square, circle, square, triangle, diamond)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c(blue,"black", "white", purple , purple, purple)[as.factor(metadat2$compartment_BCAT)])

legend("bottom",legend=c("Flow cyto control", "Bulk soil total DNA", " PCR control",  "Rhizosphere BONCAT_Active" , "Rhizosphere Inactive",  "Rhizosphere Total DNA"), 
       pch=c(15,5, 0, 1,2,5),
       cex=1.1, 
       col=c("black", "#739AFF",  "black", "#785EF0", "#785EF0",  "#785EF0"),
       bty = "n")

dev.off()









##------- log fold change -------#######
# compare active rhizosphere to active in plant #
# shows what microbes are more active in the plant than the soil

## select active taxa
ps1<- subset_samples(ps, Fraction=="BONCAT_Active" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
# 1903 taxa from non rareified data

## select active taxa in rhizosphere
ps.soil<- subset_samples(ps1,Fraction =="BONCAT_Active" & Compartment=="Rhizosphere")
otu.soil<-as.data.frame(t(as.data.frame(otu_table(ps.soil))))
head(otu.soil)
colnames(otu.soil)
# samples C10R, C1R, C2R, C5R

# select active taxa in plant
ps.plant<- subset_samples(ps1,Fraction =="BONCAT_Active" & Compartment!="Rhizosphere")
otu.plant<-as.data.frame(t(as.data.frame(otu_table(ps.plant))))
head(otu.plant)
colnames(otu.plant)
# samples C10E, C10N, C1E, C1N, C2E, C2E, C2N, C5E, C5N, C7E, C7N  
# drop sample C7 because it doesnt have a pair in the soil
otu.plant<-subset(otu.plant, select = -c(C7E.POS_S34,C7N.POS_S64))

# make df of active in rhizosphere to divid by
df<-cbind(otu.soil, otu.soil)
df<-df[sort(colnames(df))]
colnames(df)

# add 1 to everything (absent in active & absent in total = no change)

otu_log2<-log2(otu.plant/(df+.0001))

dim(otu_log2)
colnames(otu_log2)
n<-c("C10E", "C10N","C1E" , "C1N" , "C2E" , "C2N" , "C5E" , "C5N")
colnames(otu_log2)<-n
head(otu_log2)
# check distribution
otu_log2
hist(otu_log2$C10N) # most things didn't change because most taxa not present 
hist(otu_log2$C5N)
hist(otu_log2$C1E, breaks = 20)
# col dendogram doesn't work because a ton of stuff equals 0 infinit
# we re code -inf as -99
# this means a taxa wasn't present in the location sampled


otu_log2[otu_log2== "-Inf"] <- -99


############## make heatmap #1-------#################

# grab phylogeny and Read in the tree file 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree = read.tree("otu_output/tree.nwk")
# shorten phylogeny to match what is in our log 2 file
otus<-row.names(otu.raw) # all the otus
#otus we want
otus_1903<-row.names(otu_log2) 
asvs_remove<-setdiff(otus, otus_1903) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
tree.short

# order the table so it matches the 
# make otu column 
otu_log2$otus <- otus_1903

# grab correct order 
target<-tree.short$tip.label
otu_log2<-otu_log2[match(target, otu_log2$otus),]
# rm out column
otu_log2<-subset(otu_log2, select = -otus)
colnames(otu_log2)
head(otu_log2)
tree.short

# check out disctribution
#otu_log2[otu_log2=="-99"] <- NA
#summary(otu_log2)


# make matrix
     m<-as.matrix(otu_log2)
      row.names(m)
      #row.names(m)<-f # make row names family names
      # Create the matrix and get the column dendrogram for the heatmap from it.
      m<-structure(m)
      #make a dendrogram              
      col_dendro = as.dendrogram(hclust(dist(t(m))))
      
 # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
  #make a plot without phylogeny 
        #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
# focus on just taxa that are present in the plant
         
         
        
         df<-subset(otu_log2, otu_log2$C10E != -99 | otu_log2$C1E != -99  | otu_log2$C2E != -99 | otu_log2$C5E != -99 )
         # these are the ones we don't want 
         df1<-subset(otu_log2, otu_log2$C10E == -99 & otu_log2$C1E == -99  & otu_log2$C2E == -99 & otu_log2$C5E == -99)
         
         
         
         #make tree smaller
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_300<-row.names(df) 
         asvs_remove<-setdiff(otus, otus_300) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         df$otus <- otus_300
         
         # grab correct order 
         target<-tree.short$tip.label
         df<-df[match(target, df$otus),]
         # rm out column
         df<-subset(df, select = -otus)
         colnames(df)
         head(df)
         tree.short
         
         
        #make matrix      
         m<-as.matrix(df)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
  # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
 #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
         

#____________________________________________#


#### value in active compared to total for that respective location
### heat map # 2

## select total cells fraction 
  ## remove taxa from bulk soil
  ps1<- subset_samples(ps,Fraction !="Total_DNA"& Fraction!="beads" & Fraction !="ctl" ) 
  ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
  any(taxa_sums(ps1) == 0)
  ps1
  # 3173 taxa from non rareified data
         
# subset total cells fraction
    ps.Total<- subset_samples(ps1,Fraction =="Total_Cells" )
    otu.total<-as.data.frame(t(as.data.frame(otu_table(ps.Total))))
    colnames(otu.total)
    # samples in otu.total C10N, C10R, C1E, C1N, C1R, C2E, C2N, C2R, 
                          # C5E, C5R, C7E, C7N, C7R
                          #mising: C10E, C5N
                          # remove sample thats missing in active
   otu.total <- subset(otu.total, select = -C7R.SYBR_S19)
   dim(otu.total)
   colnames(otu.total)
# subset for active cells fraction    
    ps.Active<- subset_samples(ps1,Fraction =="BONCAT_Active" )
    otu.active<-as.data.frame(t(as.data.frame(otu_table(ps.Active))))
         dim(otu.active)
         colnames(otu.active)
    # missing : C7R
    # remove samples that you don't have a pair
    otu.active<-subset(otu.active, select = -c(C10E.POS_S60, C5N.POS_S63))
    dim(otu.active)
    colnames(otu.active)
 # log fold change from active to inactive 
 # add 1 to everything absent in total = no change)
         
         otu_log2<-log2(otu.active/(otu.total+.0001))
         
         dim(otu_log2)
         colnames(otu_log2)
         n<-c("C10N", "C10R" ,"C1E" , "C1N" , "C1R" , "C2E" , "C2N" , "C2R" , "C5E" , "C5R" , "C7E" , "C7N" )
         colnames(otu_log2)<-n
         head(otu_log2)
         # check distribution
         #otu_log2
         hist(otu_log2$C10N) # most things didn't change because most taxa not present 
         hist(otu_log2$C10R)
         hist(otu_log2$C1E, breaks = 20)
         # this means a taxa wasn't present in the location sampled
         # col dendogram doesn't work because a ton of stuff equals 0 infinit
         # we re code -inf as -99
         # this means a taxa wasn't present in the location sampled
         otu_log2[otu_log2== "-Inf"] <- -99
         
         
         # check out disctribution
        # otu_log2[otu_log2=="-99"] <- NA
        # summary(otu_log2)
         # change NA back to neg 99
        # otu_log2[otu_log2=="NA"] <- -99

         
         
         
########### make heatmap #2-------#################
     
         # grab phylogeny and Read in the tree file 
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
         tree = read.tree("otu_output/tree.nwk")
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_3000<-row.names(otu_log2) 
         asvs_remove<-setdiff(otus, otus_3000) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         otu_log2$otus <- otus_3000
         # grab correct order 
         target<-tree.short$tip.label
         otu_log2<-otu_log2[match(target, otu_log2$otus),]
         # rm out column
         otu_log2<-subset(otu_log2, select = -otus)
         # check names
         colnames(otu_log2)
         head(otu_log2)
         tree.short
         
         # make matrix
         m<-as.matrix(otu_log2)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
         
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))

###subset by just rhizosphere
         colnames(otu_log2)
         df<-subset(otu_log2, select = c(C10R, C1R, C2R, C5R ))
         
         
         df<-subset(df, df$C10R != -99 | df$C1R != -99  | df$C2R != -99 | df$C5R != -99 )
         # these are the ones we don't want 
         df1<-subset(df, otu_log2$C10R == -99 & otu_log2$C1R == -99  & otu_log2$C2R == -99 & otu_log2$C5R == -99)
         
         #make tree smaller
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_r<-row.names(df) 
         asvs_remove<-setdiff(otus, otus_r) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         df$otus <- otus_r
         
         # grab correct order 
         target<-tree.short$tip.label
         df<-df[match(target, df$otus),]
         # rm out column
         df<-subset(df, select = -otus)
         colnames(df)
         head(df)
         tree.short
         
         
         #make matrix      
         m<-as.matrix(df)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
  ### no phylogeny
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
         
         
### just endophyte
         colnames(otu_log2)
         df<-subset(otu_log2, select = c(C1E, C2E, C5E,  C7E))
         df<-subset(df, df$C1E != -99 | df$C2E != -99  | df$C5E != -99 | df$C7E != -99 )
         #make tree smaller
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_r<-row.names(df) 
         asvs_remove<-setdiff(otus, otus_r) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         df$otus <- otus_r
         
         # grab correct order 
         target<-tree.short$tip.label
         df<-df[match(target, df$otus),]
         # rm out column
         df<-subset(df, select = -otus)
         colnames(df)
         head(df)
         tree.short
         
         
         #make matrix      
         m<-as.matrix(df)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         ### no phylogeny
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
  ### just nodule
         colnames(otu_log2)
         df<-subset(otu_log2, select = c(C10N, C1N, C2N,  C7N))
         df<-subset(df, df$C10N != -99 | df$C1N != -99  | df$C2N != -99 | df$C7N != -99 )
         #make tree smaller
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_r<-row.names(df) 
         asvs_remove<-setdiff(otus, otus_r) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         df$otus <- otus_r
         
         # grab correct order 
         target<-tree.short$tip.label
         df<-df[match(target, df$otus),]
         # rm out column
         df<-subset(df, select = -otus)
         colnames(df)
         head(df)
         tree.short
         
         
         #make matrix      
         m<-as.matrix(df)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         ### no phylogeny
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
#######heat map #3 #######
#Are the taxa present in the plant also active in the soil?
         #active / total
         # filter for taxa that are present in the soil
         
    ### just nodule
         colnames(otu_log2)
         df<-subset(otu_log2, otu_log2$C1E != -99 | otu_log2$C2E != -99  | otu_log2$C5E != -99 | otu_log2$C7E != -99 )
         #make tree smaller
         dim(df)
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_r<-row.names(df) 
         asvs_remove<-setdiff(otus, otus_r) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         df$otus <- otus_r
         
         # grab correct order 
         target<-tree.short$tip.label
         df<-df[match(target, df$otus),]
         # rm out column
         df<-subset(df, select = -otus)
         colnames(df)
         head(df)
         tree.short
         
         
         #make matrix      
         m<-as.matrix(df)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         ### no phylogeny
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
         
#####heat map # 4######
# what is the selectivity of the plant?
# total viable in fraction X / total viable in rhizosphere
         
         
         ## select total taxa
         ps1<- subset_samples(ps, Fraction=="Total_Cells" ) 
         sample_data(ps)
         ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
         any(taxa_sums(ps1) == 0)
         ps1
         # 2473 taxa from non rareified data
##### start up <3         
         ## select active taxa in rhizosphere
         ps.soil<- subset_samples(ps1,  Compartment=="Rhizosphere")
         otu.soil<-as.data.frame(t(as.data.frame(otu_table(ps.soil))))
         head(otu.soil)
         colnames(otu.soil)
         # samples C10R, C1R, C2R, C5R, C7R
         
         # select total viable cells in plant
         ps.plant<- subset_samples(ps1,Fraction =="Total_Cells" & Compartment!="Rhizosphere")
         otu.plant<-as.data.frame(t(as.data.frame(otu_table(ps.plant))))
         head(otu.plant)
         colnames(otu.plant)
         # samples C10E, C10N, C1E, C1N, C2E, C2E, C2N, C5E, C5N, C7E, C7N  
         #otu.plant<-subset(otu.plant, select = -c(C7E.POS_S34,C7N.POS_S64))
         
         # make df of active in rhizosphere to divide by
         df<-cbind(otu.soil, otu.soil)
         df<-df[sort(colnames(df))]
         colnames(df)
         df<-subset(df, select= c(-C10R.SYBR_S20, -C5R.SYBR_S18))
         dim(df)
         dim(otu.plant)
         # add 1 to everything (absent in active & absent in total = no change)
         
         otu_log2<-log2(otu.plant/(df+.0001))
         
         dim(otu_log2)
         colnames(otu_log2)
         n<-c("C10N", "C1E","C1N", "C2E" , "C2N" , "C5E" , "C7E", "C7N")
         colnames(otu_log2)<-n
         head(otu_log2)
         # check distribution
         otu_log2
         hist(otu_log2$C10N) # most things didn't change because most taxa not present 
         hist(otu_log2$C5N)
         hist(otu_log2$C1E, breaks = 20)
         # col dendogram doesn't work because a ton of stuff equals 0 infinit
         # we re code -inf as -99
         # this means a taxa wasn't present in the location sampled
         
         
         otu_log2[otu_log2== "-Inf"] <- -11
         
         
# make heatmap
         # grab phylogeny and Read in the tree file 
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
         tree = read.tree("otu_output/tree.nwk")
         # shorten phylogeny to match what is in our log 2 file
         otus<-row.names(otu.raw) # all the otus
         #otus we want
         otus_1903<-row.names(otu_log2) 
         asvs_remove<-setdiff(otus, otus_1903) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         tree.short
         
         # order the table so it matches the 
         # make otu column 
         otu_log2$otus <- otus_1903
         
         # grab correct order 
         target<-tree.short$tip.label
         otu_log2<-otu_log2[match(target, otu_log2$otus),]
         # rm out column
         otu_log2<-subset(otu_log2, select = -otus)
         colnames(otu_log2)
         head(otu_log2)
         tree.short
         
         # check out disctribution
         #otu_log2[otu_log2=="-99"] <- NA
         #summary(otu_log2)
         
         
         # make matrix
         m<-as.matrix(otu_log2)
         row.names(m)
         #row.names(m)<-f # make row names family names
         # Create the matrix and get the column dendrogram for the heatmap from it.
         m<-structure(m)
         #make a dendrogram              
         col_dendro = as.dendrogram(hclust(dist(t(m))))
         
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(8,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
         
         
         

         
         
#____________________________________________#
         
###############------top taxa by max log fold change -------############
#select top 100 by the biggest log change 
dim(otu_log2)
otu_log2$max <- apply(abs(otu_log2), 1, max, na.rm=TRUE)

top100<-otu_log2%>% arrange(-max) %>% slice(1:100)
dim(top100)
# make otu vector
otus100<- rownames(top100)[1:100]

# remove max column
top100<-top100[,-15]
colnames(top100)
dim(top100)


# tree is so big. filter and make it smaller
# asvs we want
otus<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(otus, otus100) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
tree.short

# order the table so it matches the 
# make otu column 
top100$otus <- otus100
head(top100)

# grab correct order 
target<-tree.short$tip.label
top100<-top100[match(target, top100$otus),]
head(top100)
top100<-top100[-15] # remove otu column
colnames(top100




# get  family names
# get anesm from taxa table 
otus100
dim(taxon)
otus<-rownames(taxon)
taxon$otus <-otus
colnames(taxon)
taxa100<-subset(taxon, otus %in% otus100)
dim(taxa100)
# put taxa table in the correct order from the phylogenic tree
taxa100<-taxa100[match(target, taxa100$otus),]
# one famliy is unknown relabel that unknown
f[4]<- "unknown"
# i want f to just be a little labe across the heat map
u<-seq(1:100)
f<-paste(f, u)


# put famliy laels in 
tree.short$tip.label <-f

# make a matrix
m<-as.matrix(top100)
row.names(m)<-row.names(top100)
row.names(m)
row.names(m)<-f # make row names family names
# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(m)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
  



# And make the plot with phylogeny
#pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
windows(8,10)
heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
#dev.off()

#make a plot without phylogeny 
#make a dendrogram              
row_dendro = as.dendrogram(hclust(dist((m))))

############  make the heatmap  plot #################33
windows(8,10)
heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))




###############-----heatmap with averageing across samples ###################

#subst by treament
n<-c("C10N","C1N","C2N", "C7N") 
e<- c("C1E","C2E","C5E","C7E")
r<- c("C10R", "C1R",  "C2R",  "C5R")  
# average log fold change in each treatment
nod<- otu_log2%>% select(all_of(n))
nod<-rowSums(nod)/4
endo<- otu_log2%>% select(all_of(e))
endo<- rowSums(endo)/4
rhiz<- otu_log2%>% select(all_of(r))
rhiz<-rowSums(rhiz)/4
#combine
avg_log2<-cbind(nod, endo)%>% cbind(rhiz)

# now turn to df
avg_log2<-as.data.frame(avg_log2)


##top 100 by the biggest log change 
s<-rowSums(abs(avg_log2))
s<-sort(s, decreasing = TRUE)
top100<-names(s)[1:100]


# tree is so big. filter and make it smaller
otus<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(otus, top100) #asvs we don't want
length(otus)
length(r)
tree.endo<-drop.tip(tree, asvs_remove) # remove asvs we don't need
target<-tree.endo$tip.label

# subest dataframe to be on the top 100 taxa
# add otu col
otus<-row.names(avg_log2)
avg_log2<-cbind(avg_log2, otus)

otu_100<-subset(avg_log2, otus %in% top100)
#grab names
n<-colnames(otu_100)
r<-row.names(otu_100)
# order row to match  tree
otu_100<-otu_100[match(target, otu_100$otus),]
otu_100<-otu_100[-4] # remove otu column

# make a matrix
m<-matrix(as.numeric(unlist(otu_100)),nrow=nrow(otu_100))
row.names(m)<-r
colnames(m)<-n[-4]
#otu_log2_100

# get names of taxa
#taxon
otus<-row.names(taxon)
tax<-cbind(taxon, otus)
#subset taxa table
taxa_100<-subset(tax, otus %in% top100)
# order by otu
taxa_100<-taxa_100[match(target, taxa_100$otus),]
g<-taxa_100$Family

# rename rows
row.names(m)<-g
# rename tree tips
tree.endo$tip.label <- g

# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(m)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
# And make the plot....
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s")
svg(file="heatmap.svg",width = 8 ,height=8)
#windows(15,10)
heatmap.phylo(x = m, Rowp = tree.endo, Colp = as.phylo(as.hclust(col_dendro)))
dev.off()




######------aggregate by higher taxa level for heatmap------------########

o<-row.names(taxon) # make column to join dfs 
taxon<-cbind(taxon, o)

o<-row.names(otus)  # make column to join dfs 
otus<-cbind(otus, o)

df<-left_join(otus, taxon) #join df 
df1<-df[,c(1:42, 48)] %>%group_by(Family) %>% summarise_all(sum) #summarize by family

#make total df
total<-df1 %>% select(contains('SYBR')) #select total cells columns
total<-select(total, -C7R.SYBR_S19) #remove this sample because I don't have matching active sample for it

#make active df
active<-df1 %>% select(contains('POS')) #select total cells columns
active<-select(active, -C10E.POS_S60, -C5N.POS_S63) ## remove sample I don't have a value for inactive 

#calculate log fold change
otu_log2<-log2(active+1/(total+1))

#label columns and rows
n<-c("C10N", "C10R" ,"C1E" , "C1N" , "C1R" , "C2E" , "C2N" , "C2R" , "C5E" , "C5R" , "C7E" , "C7N" )
colnames(otu_log2)<-n
row.names(otu_log2)<-df1$Family 
head(otu_log2)


######avg by treament############
n<-c("C10N","C1N","C2N", "C7N") 
e<- c("C1E","C2E","C5E","C7E")
r<- c("C10R", "C1R",  "C2R",  "C5R")  
# average log fold change in each treatment
nod<- otu_log2%>% select(all_of(n))
nod<-rowSums(nod)/4
endo<- otu_log2%>% select(all_of(e))
endo<- rowSums(endo)/4
rhiz<- otu_log2%>% select(all_of(r))
rhiz<-rowSums(rhiz)/4
#combine
avg_log2<-cbind(nod, endo)%>% cbind(rhiz)

# now turn to df
avg_log2<-as.data.frame(avg_log2)


#########top 100 by the biggest log change ############
s<-rowSums(abs(avg_log2))
s<-sort(s, decreasing = TRUE)
top100<-names(s)[1:100]


## Read in the tree file made by jenn
## rip this tre is by otus
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree = read.tree("otu_output/tree.nwk")
print(tree)

###### lets do something jenky ###
# grab a representative otu for each family.

df # df of outs and taxa
df1<-df%>% filter(Family %in% top100) # select top familys
sums<-rowSums(df1[,c(1:42)]) # sum across samples
df1<-data.frame(df1$o, df1$Family, sums) # make new df
df1<-(df1[order(df1$df1.Family, -df1$sums),]) # order df
# grab most abundance otu for each family
df1.o=1
df1.Family=1
sums=1
newdf <- data.frame(df1.o, df1.Family, sums)
df1<-filter(df1, df1$df1.Family != "") # drop famliy that is missing a label
top100<-unique(df1$df1.Family)
  
for (i in top100) {
  print(i)
  d<- filter(df1, df1$df1.Family == i) # select df
  d<-(d[order(-d$sums),]) # order df
  d<-d[1,] # grab top taxa
  print(d)
  newdf<-rbind(d, newdf)
}

rep.otus<-newdf$df1.o


# tree is so big. filter and make it smaller
otus<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(otus, rep.otus) #asvs we don't want
length(otus)
asvs_remove
tree.endo<-drop.tip(tree, asvs_remove) # remove asvs we don't need
target<-tree.endo$tip.label
print(tree.endo)
### change names in tre to be family names
 # order otus in df that matches otus and familys
newdf1<-newdf[match(target, newdf$df1.o),]
newdf1
# put family names on tree
#taxon
g<-newdf1$df1.Family
# rename tree tips
tree.endo$tip.label <- g

# rename rows
row.names(m)<-g





# subest avg_log2 dataframe to be on the top 100 taxa

# add famliy col

fam<-row.names(avg_log2)
avg_log2<-cbind(avg_log2, fam)

otu_100<-subset(avg_log2, fam %in% top100)
#grab names
dim(otu_100)
n<-colnames(otu_100[1:3])
r<-row.names(otu_100)

# order row to match  tree
otu_100<-otu_100[match(target, otu_100$otus),]
otu_100<-otu_100[-4] # remove otu column
otu_100<-otu_100[match(newdf1$df1.Family, otu_100$fam),]

# make a matrix
m<-matrix(as.numeric(unlist(otu_100[1:3])),nrow=nrow(otu_100))
row.names(m)<-r
colnames(m)<-n[-4]
#otu_log2_100


# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(m)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
# And make the plot....
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s")
svg(file="heatmap.svg",width = 8 ,height=8)
#windows(15,10)
heatmap.phylo(x = m, Rowp = tree.endo, Colp = as.phylo(as.hclust(col_dendro)))
dev.off()




#################------ fold change total to active---------------####
#### remove taxa from bulk soil
ps1<- subset_samples(ps,Fraction !="Total_DNA"& Fraction!="beads" & Fraction !="ctl" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
# 6189 taxa

#Who are the most abundant taxa minus the total DNA
#who are most abundant taxa
topN = 200
most_abundant_taxa = sort(taxa_sums(ps1), TRUE)[1:topN]
print(most_abundant_taxa)
Top_tax = prune_taxa(names(most_abundant_taxa), ps1)
length(get_taxa_unique(Top_tax, "Class"))
length(get_taxa_unique(Top_tax, "Phyla"))
length(get_taxa_unique(Top_tax, "Genus"))
length(get_taxa_unique(Top_tax, "Family"))
tax_table(Top_tax)
otu_table(Top_tax)



# subset total cell and active fraction
ps.Total<- subset_samples(ps1,Fraction =="Total_Cells" )
otu_total<-as.data.frame(t(as.data.frame(otu_table(ps.Total))))
# boncat pos table doesn't have C7Rpos, so I'll remove that from the table

head(otu_total)


otu_total<-select(otu_total, -C7R.SYBR_S19)
# length #7185
#       C10N.SYBR_S26 C10R.SYBR_S20 C1E.SYBR_S21 C1N.SYBR_S13 
# taxa1  139581          1430        62828       126909  
# taxa2  19537           179         8553        17541
# taxa3  78                0         4336        22797

ps.Active<- subset_samples(ps1,Fraction =="BONCAT_Active" )
otu_active<-as.data.frame(t(as.data.frame(otu_table(ps.Active))))
## missing some sample from sybr for C10 E and C5 N
## so i think I'll just remove those ones and won't have a value for inactive for those samples
otu_active<-select(otu_active, -C10E.POS_S60, -C5N.POS_S63)

#head(otu_active)
#sum<-rowSums(otu_active)
# length 7185
# taxa1 10000
# taxa2 0
# taxa3 0 


##------- log fold change from active to inactive ----------

#what do you do for things that are not present in total?
# can we just add 1 to everything?

otu_log2<-log2(otu_active+1/(otu_total+1))


n<-c("C10N", "C10R" ,"C1E" , "C1N" , "C1R" , "C2E" , "C2N" , "C2R" , "C5E" , "C5R" , "C7E" , "C7N" )
colnames(otu_log2)<-n
head(otu_log2)
#zero in numerator = not present in active = -inf

# check distribution
otu_log2
hist(otu_log2$C10N)
hist(otu_log2$C10R)
hist(otu_log2$C1E, breaks = 20)

#####------wrangling log2 data -------#####

#transpose
t_log2_otu<- t(otu_log2)

#import metadata for summarising 
metadat_l<- read.delim("16s/metadata_log2.txt", sep="\t", header = T, check.names=FALSE)
y<-colnames(otu_log2)
rownames(metadat_l) <- y
metadat_l<-as.data.frame(metadat_l)
#insert metadata to data frame
t_log2_otu<-cbind(metadat_l, t_log2)
#sumarise by rep
t_log2_otu_short<- t_log2_otu %>% group_by(Fraction) %>% summarise_if(is.numeric,median) %>% select(-REP)

# filter dataset for top asvs
# I selectws the top taxa by relative abundance above

otu_table(Top_tax)
df<-as.data.frame(tax_table(Top_tax))
dim(df)
#I'm curious who these asvs are so I'm going to add there names in
asvs<-row.names(df)
df<-mutate(df, asvs=asvs)
asvs<-row.names(otu_log2)
otu_log2<-mutate(otu_log2, asvs=asvs)

#filter the log2 dataset by the top taxa
otu_log2_200<-left_join(df, otu_log2, by= "asvs")

as.matrix(otu_log2_200[,-(1:8)])


## Read in the tree file made by jenn
tree = read.tree("16s/tree.nwk")
# my tree just has ASV names think

et the column dendrogram for the heatmap from it.
# got remove all those extract columns and 

asvs200 =   (otu_log2_200[, 'asvs'])
m = as.matrix(otu_log2_200[,-(1:8)])
row.names(m)<-asvs200 
m<-structure(m)

#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))

# lol tree is so big ? maybe filter and make it smaller
# we need a list of all the asvs we don't want

asvs200 # asvs we want

dim(otus.raw)

asvs<-row.names(otus.raw)

asvs_remove<-setdiff(asvs, asvs200)


tree<-drop.tip(tree, asvs_remove)

# And )make the plot....
#pdf(file="../figures/Fig2_heatmap24wk_Phylo.pdf",width = 5,height=8, useDingbats=FALSE)
heatmap.phylo(x = m, Rowp = tree, Colp = as.phylo(as.hclust(col_dendro)))
#dev.off()





# select top taxa but select the taxa with the biggest change?



#import it phyloseq
Workshop_OTU <- otu_table(as.matrix(otus.phyloseq), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat_l)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
psl <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

#test it worked
sample_names(psl)
print(psl)
# 7185 taxa because we didn't include anything that wasn't rhizo, nodule or end












#---- calculating inactive fraction-------
#Inactive<- total- active
otu_inactive<-otu_total-otu_active
# all the taxa that are higher in the active fraction than the inactive fraction (negative )are == to zero 
otu_inactive[otu_inactive<0]<-0
otu_inactive
#length 7185
#       C10N.SYBR_S26 C10R.SYBR_S20 C1E.SYBR_S21 C1N.SYBR_S13 
# taxa1  0              47          0            0  
# taxa2  0              0           0            0
# taxa3  0              22          0            0
# change colnames
n<-c("C10N.inactive", "C10R.inactive", "C1E.inactive"  ,"C1N.inactive" , "C1R.inactive", "C2E.inactive" , "C2N.inactive"  ,"C2R.inactive" , "C5E.inactive", 
"C5R.inactive",  "C7E.inactive",  "C7N.inactive" )
colnames(otu_inactive)<- n
# append the original dataset by adding the inactive taxa
#make an asvs col to join by
y<-row.names(otus)
otus<-mutate(otus, asvs= y)

y<-row.names(otu_inactive)
otu_inactive<-mutate(otu_inactive, asvs= y)

# join data frames
otus<-full_join(otu_inactive, otus)
# as row names from the asvs bit
y<-otus$asvs
row.names(otus) <- y
otus<-select(otus, -asvs)
otus
# make all the asvs that are NA in the inactive fraction zeros
otus[is.na(otus)] <- 0


# append metadata
setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")
metadat <- read.delim("16s/metadata_w_inactive.txt", sep="\t", header = T, check.names=FALSE)
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule', 'ctl'='ctl'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))

#to make coloring things easier I'm gong to added a combined fractionXboncat column no sure if i need this
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))
metadat<-select(metadat, -BONCAT)

## order metadata
metadat<-metadat[order(metadat$SampleID),]
metadat
## order otu table
otus<-otus[, order(colnames(otus))]

## Convert OTU numbers to percentages ##
otus.perc<-otus/rowSums(otus)*100


##### Phyloseq and manipulation by taxonomy ####

#check n taxa
rank_names(ps)

#what phyla are here?
rank_names(ps)
get_taxa_unique(ps, "Phyla")
#how many of each taxa?
taxa_sums(ps)

# Are some taxa rare?
OTUabundance = as.numeric(taxa_sums(ps))
ASV = taxa_names(ps)
d1<-as.data.frame(OTUabundance)
d1<-cbind(d1, ASV)
windows()
ggplot(d1, aes(OTUabundance)) + 
  geom_histogram() +
  ggtitle("Histogram of Total Counts") + 
  xlim(0, 20000) + ylim (0,50) + theme_bw()
# most taxa are rare, but some are really abundant

# who are the abundant taxa?
topN = 50
most_abundant_taxa = sort(taxa_sums(ps), TRUE)[1:topN]
print(most_abundant_taxa)
GP20 = prune_taxa(names(most_abundant_taxa), ps)
length(get_taxa_unique(GP20, "Class"))
print(get_taxa_unique(GP20, "Phyla"))
print(get_taxa_unique(GP20, "Family"))
abundant_tax<-tax_table(GP20)

# who is in the negative PCR ctl?
#subset
ctl<-subset_samples(ps, SampleID=="CTL_S66")
sample_variables(ps)
# rm taxa that aren't there
any(taxa_sums(ctl) == 0)
ctl<-prune_taxa(taxa_sums(ctl) > 0, ctl)
any(taxa_sums(ctl) == 0)
ntaxa(ctl)
#who are the most abundant taxa?
topN = 20
most_abundant_taxa = sort(taxa_sums(ctl), TRUE)[1:topN]
print(most_abundant_taxa)
ctl20 = prune_taxa(names(most_abundant_taxa), ctl)
length(get_taxa_unique(ctl20, "Family"))
print(get_taxa_unique(ctl20, "Phyla"))
get_taxa_unique(ctl20, "Genus")
get_taxa_unique(ctl20, "Family")

#who is in the flow cytometer ctl?
#subset
ctl<-subset_samples(ps, SampleID=="BEADS_S67")
sample_variables(ps)
# rm taxa that aren't there
any(taxa_sums(ctl) == 0)
ctl<-prune_taxa(taxa_sums(ctl) > 0, ctl)
any(taxa_sums(ctl) == 0)
ntaxa(ctl)
#who are the most abundant taxa?

topN = 20
most_abundant_taxa = sort(taxa_sums(ctl), TRUE)[1:topN]
print(most_abundant_taxa)
ctl20 = prune_taxa(names(most_abundant_taxa), ctl)
length(get_taxa_unique(ctl20, "Family"))
print(get_taxa_unique(ctl20, "Phyla"))
get_taxa_unique(ctl20, "Genus")
get_taxa_unique(ctl20, "Family")

#who is in the nodule?
#subset
nod<-subset_samples(ps, Compartment=="Nodule")

# rm taxa that aren't there
any(taxa_sums(nod) == 0)
nod<-prune_taxa(taxa_sums(nod) > 0, nod)
any(taxa_sums(nod) == 0)
ntaxa(nod)

#who are most abundant taxa
topN = 20
most_abundant_taxa = sort(taxa_sums(nod), TRUE)[1:topN]
print(most_abundant_taxa)
nod20 = prune_taxa(names(most_abundant_taxa), nod)
length(get_taxa_unique(nod20, "Class"))
print(get_taxa_unique(nod20, "Phyla"))
get_taxa_unique(nod20, "Genus")
get_taxa_unique(nod20, "Family")
get_taxa_unique(nod20, "Class")
tax_table(nod20)


######3.------  import percent abundance into phyloseq for figure ----- #####

otus.phyloseq<- t(otus.perc)

#import it phyloseq
Workshop_OTU <- otu_table(as.matrix(otus.phyloseq), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps_perc <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

t(otu_table(ps_perc))

#test it worked
sample_names(ps_perc)
print(ps_perc)
# 14593 taxa

topN = 50
most_abundant_taxa = sort(taxa_sums(ps_perc), TRUE)[1:topN]
print(most_abundant_taxa)
ps_20 = prune_taxa(names(most_abundant_taxa), ps_perc)
length(get_taxa_unique(ps_perc, "Class"))
print(get_taxa_unique(ps_20, "Phyla"))
get_taxa_unique(ps_perc, "Genus")
get_taxa_unique(ps_perc, "Family")
get_taxa_unique(ps_perc, "Class")
tax_table(ps_20)

# top phyla

#" Cyanobacteria"     " Acidobacteriota"   " Armatimonadota"    " Nitrospirota"      " Proteobacteria"    " Gemmatimonadota"   " Chloroflexi"      
#" Bacteroidota"      " Actinobacteriota"  " Planctomycetota"   " Verrucomicrobiota" " Myxococcota"       " Bdellovibrionota" 

#### 100% plots

Cyanobacteria <- subset_taxa(ps, Phyla ==  " Cyanobacteria"  )
cyano.sum<-rowSums(otu_table(Cyanobacteria))
Acidobacteriota <- subset_taxa(ps, Phyla  == " Acidobacteriota" )
Acido.sum<-rowSums(otu_table(Acidobacteriota))
Armatimonadota <- subset_taxa(ps, Phyla  == " Armatimonadota"  )
Arma.sum<-rowSums(otu_table(Armatimonadota))
Nitrospirota <- subset_taxa(ps, Phyla  == " Nitrospirota"   )
Nitro.sum<-rowSums(otu_table(Nitrospirota))
Proteobacteria  <- subset_taxa(ps, Phyla == " Proteobacteria"  )
Proteo.sum<-rowSums(otu_table(Proteobacteria))
Gemmatimonadota <- subset_taxa(ps, Phyla  ==  " Gemmatimonadota"  )
Gemma.sum<-rowSums(otu_table(Gemmatimonadota))
Chloroflexi <- subset_taxa(ps, Phyla  ==    " Chloroflexi" )
Chloro.sum<-rowSums(otu_table(Chloroflexi))
Bacteroidota  <- subset_taxa(ps, Phyla  ==    " Bacteroidota"  )
Bactero.sum<-rowSums(otu_table(Bacteroidota))
Actinobacteriota  <- subset_taxa(ps, Phyla  ==     " Actinobacteriota"   )
Actino.sum<-rowSums(otu_table(Actinobacteriota ))
Planctomycetota <- subset_taxa(ps, Phyla  ==     " Planctomycetota" )
Planctomycetota.sum<-rowSums(otu_table(Planctomycetota ))
Verrucomicrobiota  <- subset_taxa(ps, Phyla  ==     " Verrucomicrobiota"   )
Verrucomicrobiota.sum<-rowSums(otu_table(Verrucomicrobiota ))
Myxococcota <- subset_taxa(ps, Phyla  ==     " Myxococcota"   )
Myxococcota.sum<-rowSums(otu_table(Myxococcota ))
Bdellovibrionota <- subset_taxa(ps, Phyla  ==      " Bdellovibrionota"   )
Bdellovibrionota.sum<-rowSums(otu_table(Bdellovibrionota ))


#other<-100-(Bdellovibrionota.sum+ Myxococcota.sum+ Verrucomicrobiota.sum + Planctomycetota.sum+ Actino.sum + Bactero.sum +
#            cyano.sum + Acido.sum + Arma.sum + Nitro.sum + Proteo.sum + Gemma.sum + Chloro.sum)


phyl.mat<-cbind(Bdellovibrionota.sum, Myxococcota.sum, Verrucomicrobiota.sum , Planctomycetota.sum, Actino.sum , Bactero.sum ,
                  cyano.sum , Acido.sum , Arma.sum , Nitro.sum , Proteo.sum , Gemma.sum , Chloro.sum) 
  
phyl.perc<-(phyl.mat/41610)*100  
print(phyl.perc)
other <- 100-rowSums(phyl.perc)
phyl.perc<- cbind(phyl.perc, other)
phyl.perc<-as.data.frame(phyl.perc)
phyl.perc<-cbind(metadat, phyl.perc)


# data wrangling  
phy.df<-gather(phyl.perc, "taxa", value, 7:20 )

head(phy.df)

df<-phy.df %>% filter(Fraction!="Inactive")

df<-df%>%
  group_by(taxa, compartment_BCAT) %>%
  summarise(mean = mean(value))

unique(df$compartment_BCAT)

df$compartment_BCAT<-factor(df$compartment_BCAT, levels = c("beadsbeads", 
                    "ctlctl", "Bulk_SoilTotal_DNA",  "RhizosphereTotal_DNA",
                    "RhizosphereTotal_Cells", "RhizosphereBONCAT_Active" ,"RootsTotal_Cells"  
                    , "RootsBONCAT_Active", "NoduleTotal_Cells"  ,  "NoduleBONCAT_Active"   ))

df<-filter(df, compartment_BCAT!="beadsbeads" & compartment_BCAT!="ctlctl")
# Stacked + percent
svg(file="top_taxa.svg",width = 6, height=5 )
install.packages("viridis")  # Install
library("viridis")
windows()

df%>% 
ggplot(aes(fill=taxa, y=mean, x=compartment_BCAT)) + 
  geom_bar(position="fill", stat= "identity")+
  #scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Top phyla") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Compartment")+
  ylab("relative abundance %")
  
dev.off()






######------- 100% plot for just active taxa---------

#maybe try subsetting by sample?
ps_active<-subset_samples(ps, Fraction !=  "Total_Cells")


metadat_act<-subset(metadat, Fraction!="Total_Cells")

# 100% plots 

rhizobiaceae <- subset_taxa(ps_active, Family ==  " Rhizobiaceae"  )
Rhizobiaceae<-rowSums(otu_table(rhizobiaceae))
psuedo <- subset_taxa(ps_active, Family == " Pseudomonadaceae" )
Pseudomonadaceae<-rowSums(otu_table(psuedo))
staph <- subset_taxa(ps_active, Family == " Staphylococcaceae")
Staphylococcaceae<-rowSums(otu_table(staph))
burk <- subset_taxa(ps_active, Family == " Burkholderiaceae"  )
Burkholderiaceae<-rowSums(otu_table(burk))
sphing <- subset_taxa(ps_active, Family== " Sphingomonadaceae" )
Sphingomonadaceae<-rowSums(otu_table(sphing))
micro <- subset_taxa(ps_active, Family ==  " Micrococcaceae"   )
Micrococcaceae<-rowSums(otu_table(micro))

other<-100-(Rhizobiaceae+Pseudomonadaceae +Staphylococcaceae +Burkholderiaceae + Sphingomonadaceae + Micrococcaceae)




phyl.mat<-cbind(Rhizobiaceae, Pseudomonadaceae, Staphylococcaceae, Burkholderiaceae,  Sphingomonadaceae , Micrococcaceae, other)
print(phyl.mat)
plyl.mat<-as.data.frame(phyl.mat)
phyl.mat<-aggregate(phyl.mat~metadat_act$compartment_BCAT,FUN=mean)

# data wrangling  
phy.mat<-gather(phyl.mat, "taxa", value, 2:8)
colnames(phy.mat)[1]<-"compartment_BCAT"

phy.df<-phy.mat %>% group_by(compartment_BCAT, taxa) %>%
  summarise(mean = mean(value), n = n())

phy.df$compartment_BCAT<-factor(phy.df$compartment_BCAT, levels = c("NActl", "Bulk_SoilTotal_DNA",  "RhizosphereTotal_DNA","RhizosphereBONCAT_Active" ,"RootsBONCAT_Active",  "NoduleBONCAT_Active"   ))



# Stacked + percent

svg(file="figures/16s/top_active_taxa.svg",width = 8, height=6 )

#windows()
ggplot(phy.df, aes(fill=taxa, y=mean, x=compartment_BCAT)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  ggtitle("Top Families") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  xlab("Compartment")+
  ylab("relative abundance %")

dev.off()

#------- 3. Run PCoA analysis of entire dataset on Otus------


# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps), method = "bray")

# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(54-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

#------------plant verse soil---------

# subset
p<-otus.p[which(metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive"),]

# subset metadata
metadat_in<-metadat%>% filter(Fraction !="Total_Cells" & metadat$Fraction != "Inactive" ) 
as.factor(metadat_in$Compartment)
as.factor(metadat_in$Fraction)
as.factor(metadat_in$compartment_BCAT)
levels(as.factor(metadat_in$compartment_BCAT))
levels(as.factor(metadat_in$Fraction))

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/Pcoa_plantvsoil_raw.svg",width = 7, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(p, col=c("black"),
       pch=c(22,21,0, 23 )[as.factor(metadat_in$Fraction)],
       lwd=1,cex=2,
       bg=c("black", "#739AFF", "white", "#DC267F", "#785EF0", "#785EF0" , "#FFB000")[as.factor(metadat_in$compartment_BCAT)])


legend("top",legend=c("beads", "Bulk_SoilTotal_DNA" ,      "neg control"   ,  "NoduleBONCAT_Active"    ,  "NoduleInactive"     ,      "RhizosphereBONCAT_Active",
                       "RhizosphereInactive","RhizosphereTotal_DNA" ,    "RootsBONCAT_Active"   ,    "RootsInactive"),
       pch=c(15,5,0,1,2 , 1, 2 , 5, 1, 2),
       col= c("black", "#739AFF", "black", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000"),
       bty = "n",
       inset = c(.05, 0))

dev.off()

#####-------active verse inactive---------------

# subset
p<-otus.p[which(metadat$Fraction != "Total_Cells"),]


# subset metadata
metadat_in<-metadat%>% filter(Fraction !="Total_Cells") 

metadat%>% filter(Fraction == "Inactive")
as.factor(metadat_in$Compartment)
as.factor(metadat_in$Fraction)
as.factor(metadat_in$compartment_BCAT)
levels(as.factor(metadat_in$compartment_BCAT))
levels(as.factor(metadat_in$Fraction))

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/pcoa/pcoa_all_raw.svg",width = 7, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(p, col=c("black"),
       pch=c(22,21,0,24, 23 )[as.factor(metadat_in$Fraction)],
       lwd=1,cex=2,
       bg=c("black", "#739AFF", "white", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000")[as.factor(metadat_in$compartment_BCAT)])


legend("bottom",legend=c("beads", "Bulk_SoilTotal_DNA" ,      "neg control"   ,  "NoduleBONCAT_Active"    ,  "NoduleInactive"     ,      "RhizosphereBONCAT_Active",
                      "RhizosphereInactive","RhizosphereTotal_DNA" ,    "RootsBONCAT_Active"   ,    "RootsInactive"),
       pch=c(15,5,0,1,2 , 1, 2 , 5, 1, 2),
       col= c("black", "#739AFF", "black", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000"),
       bty = "n",
       inset = c(.05, 0))

dev.off()



##########-------run a new pcoa on just soil <3--------###########
######forget this pc3 and pc 4 stufff for rn
ps
ps2<-subset_samples(ps, Compartment !=  "Nodule" & Compartment != "Roots")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")

# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(28-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots") 

as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/soil_raw.svg",width = 6, height=6 )
windows(title="PCoA on asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(22,21,0,24, 25, 23)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c("black","#739AFF", "white", "#785EF0", "#785EF0",  "#785EF0", "#785EF0" )[as.factor(metadat2$compartment_BCAT)])

# bulk soil = blue
# rhizzo = purple
# ctl = white
# bead = black
# active = circle
# inactive = trangle
legend("topleft",legend=c("Flow cyto control", "Bulk soil total DNA", " PCR control",  "Rhizosphere BONCAT_Active" , "Rhizosphere Inactive",  "Rhizosphere Total DNA"), 
       pch=c(15,5, 0, 1,2,5),
       cex=1.1, 
       col=c("black", "#739AFF",  "black", "#785EF0", "#785EF0",  "#785EF0"),
       bty = "n")

dev.off()


#########-------------- permanova----------------#########



otu.p.t <- t(otus.perc)
nrow(otu.p.t)
otu.perm<- adonis2(otu.p.t~ Compartment*Fraction, data = metadat, permutations = 999, method="bray")

otu.perm
# Fraction         4   8.3386 0.59802 18.4333  0.001 ***
#  BONCAT           2   1.2988 0.09315  5.7422  0.001 ***
#  Fraction:BONCAT  2   0.5742 0.04118  2.5387  0.126   

#analysis of similarities
otu.ano<- anosim(otu.p.t, grouping =  metadat$Compartment, permutations = 999)
summary(otu.ano)

#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=metadat$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

dispersion <- betadisper(otus.bray, group=metadat$BONCAT)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Groups     3 0.16064 0.053547 1.148    999  0.363
# homogenous varience :)

dispersion <- betadisper(otus.bray, group=metadat$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Groups     8 0.69417 0.086771 26.697    999  0.001 ***
# non homogenous variance :(


##### subset data by fraction 
#rhizo active ver inactive
test<-otu.p.t[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Total_DNA"),]
metadat_y<-metadat[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Total_DNA") ,]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm
#rhizo active ver total DNA
test<-otu.p.t[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive"),]
metadat_y<-metadat[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive") ,]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm


#nod active verse inactive
test<-otu.p.t[which(metadat$Compartment == "Nodule" &  metadat$Fraction != "Total_Cells" ),]
metadat_y<-metadat[which(metadat$Compartment == "Nodule" & metadat$Fraction != "Total_Cells" ),]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm

#endo active verse inactive
test<-otu.p.t[which(metadat$Compartment== "Roots" & metadat$Fraction != "Total_Cells" ),]
metadat_y<-metadat[which(metadat$Compartment == "Roots" & metadat$Fraction != "Total_Cells" ),]
otu.perm<- adonis2(test~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm


#active roots verse nodules
test<-otu.p.t[which(metadat$Compartment != "Rhizosphere" &  metadat$Compartment != "Bulk_Soil" & metadat$Fraction == "BONCAT_Active"),]
metadat_y<-metadat[which(metadat$Compartment != "Rhizosphere" & metadat$Compartment != "Bulk_Soil" &   metadat$Fraction == "BONCAT_Active") ,]
otu.perm<- adonis2(test~ Compartment, data = metadat_y, permutations = 999, method="bray")
otu.perm

#active roots verse rhizo
test<-otu.p.t[which(metadat$Compartment != "Rhizosphere" &  metadat$Compartment != "Bulk_Soil" & metadat$Fraction == "BONCAT_Active"),]
metadat_y<-metadat[which(metadat$Compartment != "Rhizosphere" & metadat$Compartment != "Bulk_Soil" &   metadat$Fraction == "BONCAT_Active") ,]
otu.perm<- adonis2(test~ Compartment, data = metadat_y, permutations = 999, method="bray")
otu.perm

