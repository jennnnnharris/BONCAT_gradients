# 16s analyis

### 1. Initial Setup ###

### Clear workspace ###

rm(list=ls())

## if trouble loading phyloseq, see: http://joey711.github.io/phyloseq/install

#source("http://bioconductor.org/biocLite.R")
#biocLite("Heatplus")

#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("phyloseq")
BiocManager::install("Heatplus")


### Load required libraries ###

library(gplots)
library(ggplot2)
library(vegan)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(data.table)
library(dplyr)
library(phyloseq)
library(DT)
library(Heatplus)
library(viridis)
library(hrbrthemes)
library(ade4)

## Set the working directory; modify to your own ###
setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")

### Import Data ###
taxon <- read.table("16s/taxonomy.txt", sep="\t", header=T, row.names=1)
otus <- read.table("16s/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("16s/metadata.txt", sep="\t", header = T, check.names=FALSE)

## Transpose OTU table ##
otus.t <- t(otus)

## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order otu table
otus.t<-otus.t[order(row.names(otus.t)),]

## Determine minimum available reads per sample ##
min(rowSums(otus.t))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
otus.r<-rrarefy(otus.t, 41610)

## Convert OTU numbers to percentages ##
otus.perc<-otus.r/rowSums(otus.r)*100

#set colors
mycols=c("grey27","grey","#9795ff","#4406e2", "#ffb4f6","#e20a8f", "#695e00", "#6eda0a", "#0a7416" )

# rhizo, endo, nod, 
#green , grey, purple , grey, pink, grey
mycols3= c("#3aaf04", "grey","#4406e2", "grey" , "#ff50c8", "grey")

# recode some names so they are easier to understand

#Fraction = BONCAT active, Total cells
#Compartments = 
metadat$Fraction
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat$Compartment<-factor(metadat$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active' ))

#### rarefaction curve


data(dune)
min(rowSums(t(dune)))

head(dune)
specaccum

sp<-specaccum(dune, method = "random", permutations = 100,
          conditioned =TRUE, gamma = "jack1")

sp$sites
plot(sp$sites, sp$richness)
#yay seemed to work with test data set

row.names(otus.t)

sp<-specaccum(otus.t, method = "random", permutations = 100,
              conditioned =TRUE, gamma = "jack1")
plot(sp$sites, sp$richness, xlab="sample", ylab = "# ASVs")


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
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$BONCAT))

## Plot ordination with factors coloured and shaped as you like # 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

svg(file="figures/16s/Pcoa.svg",width = 6, height=6 )
#windows(title="PCoA on OTUs - Bray Curtis")
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, col=c("black"),
       pch=c(21,21,22,23,24)[as.factor(metadat$Compartment)],
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
otu.perm<- adonis2(otus.perc~ Compartment*BONCAT, data = metadat, permutations = 999, method="bray")
otu.perm
# Fraction         4   8.3386 0.59802 18.4333  0.001 ***
#  BONCAT           2   1.2988 0.09315  5.7422  0.001 ***
#  Fraction:BONCAT  2   0.5742 0.04118  2.5387  0.126   

#analysis of similarities
otu.ano<- anosim(otus.perc, grouping =  metadat$Compartment, permutations = 999)
summary(otu.ano)

#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=metadat$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Groups     4 0.61336 0.153340 40.69    999  0.001 ***
# non homogeneous varience :( 

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
#rhizo
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

#endo
test<-otus.perc[which(metadat$Fraction == "Endo"),]
metadat_t<-metadat[which(metadat$Fraction == "Endo"),]
otu.perm<- adonis2(test~ BONCAT, data = metadat_t, permutations = 999, method="bray")
otu.perm
#BONCAT    1  0.19218 0.58686 9.9433  0.014 *

######------make phyloseq object-------#####

otus.phyloseq<- t(otus.t)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otus)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_OTU <- otu_table(as.matrix(otus.t), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

#test it worked
sample_names(ps)
print(ps)

#rarafaction curve in phyloseq

?rarecurve
  
rarecurve((otu_table(ps)), step=50, cex=0.5)


#####2. Calculate diversity#######

# diversity 
rich<-estimate_richness(ps, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
p <- plot_richness(ps, "Fraction", measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
p <- p + geom_boxplot(aes(fill = "Fraction")) + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578"))
print(p)

# set wd for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

# diversity of active microbes in each fraction
rich<-cbind(rich, metadat)
rich<-as.data.frame(rich)
colnames(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))


svg(file="figures/16s/alpha_diversity.svg",width = 5, height=4 )
rich %>%
  filter(BONCAT!="DNA", Fraction!="ctl")%>%
  ggplot(aes(x=Compartment, y=Observed, fill=Fraction,))+
  geom_boxplot() +
  scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  ylab("Number of ASVs")+
  theme_bw()
dev.off()


#phylogenic diversity
#https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#alpha-diversity-calculations


#### 3. Phyloseq and manipulation by taxonomy ####

# get rid of taxa that aren; tin any samples
Workshop.16S<-prune_taxa(taxa_sums(Workshop.16S) > 0, Workshop.16S)
any(taxa_sums(Workshop.16S) == 0)

#check n taxa
rank_names(Workshop.16S)

# remove chloroplast DNA
Workshop.16S<-subset_taxa(Workshop.16S, Class!="Chloroplast")

# We make a data table with information on the OTUs
ps1.dt.taxa = data.table(tax_table(Workshop.16S),OTUabundance = taxa_sums(Workshop.16S),OTU = taxa_names(Workshop.16S))
ps1.dt.tax.plot <- ggplot(ps1.dt.taxa, aes(OTUabundance)) + geom_histogram() + ggtitle("Histogram of OTU (unique sequence) counts") + theme_bw()
print(ps1.dt.tax.plot)
ggplot(ps1.dt.taxa, aes(OTUabundance)) + 
  geom_histogram() +
  ggtitle("Histogram of Total Counts") + 
  xlim(0, 1000) + ylim (0,50) + theme_bw()


plot_bar(Workshop.16S, "Fraction", "Abundance", "Phyla", title=title)


#interacting with phyloseq object
sample_variables(Workshop.16S)
length(sample_variables(Workshop.16S))

#what phyla are here?
rank_names(Workshop.16S)
get_taxa_unique(Workshop.16S, "Phyla")
#how many of each taxa?
taxa_sums(Workshop.16S)

#select most abundant taxa
topN = 100
most_abundant_taxa = sort(taxa_sums(Workshop.16S), TRUE)[1:topN]
print(most_abundant_taxa)
GP20 = prune_taxa(names(most_abundant_taxa), Workshop.16S)
length(get_taxa_unique(GP20, "Phyla"))
print(get_taxa_unique(GP20, "Phyla"))

#lookign at ctl
windows()
plot_bar(Workshop.16S, fill = "Phyla",facet_grid = "Fraction~.")
print(get_taxa_unique(GP20, "Family"))
sample_variables(Workshop.16S)

#subet
ctl<-subset_samples(Workshop.16S, Fraction=="ctl")

# rm taxa that aren;t there
any(taxa_sums(ctl) == 0)
ctl<-prune_taxa(taxa_sums(ctl) > 0, ctl)

ntaxa(ctl)
get_taxa_unique(ctl, "Phyla")


#### 100% plots ####

#I think this is transposed?

Actino <- subset_taxa(Workshop.16S, Phyla = "Actinobacteriota")
Actino.sum<-rowSums(otu_table(Actino))
Proteo <- subset_taxa(Workshop.16S, Phyla = "Proteobacteria")
Proteo.sum<-rowSums(otu_table(Proteo))
Acid <- subset_taxa(Workshop.16S, Phyla = "Acidobacteria")
Acid.sum<-rowSums(otu_table(Acid))
cyano <- subset_taxa(Workshop.16S, Phyla = "Cyanobacteria")
cyano.sum<-rowSums(otu_table(cyano))
Firmi <- subset_taxa(Workshop.16S, Phyla = "Firmicutes")
Firmi.sum<-rowSums(otu_table(Firmi))
Bact <- subset_taxa(Workshop.16S, Phyla = " Bacteroidota" )
Bact.sum<-rowSums(otu_table(Bact))


Other<-100-(Actino.sum+Acid.sum+Firmi.sum+Bact.sum+Proteo.sum+cyano.sum)

phyl.mat<-cbind(Actino.sum,Proteo.sum,Acid.sum,Firmi.sum,Bact.sum,cyano.sum,Other)
print(phyl.mat)
phyl.ag<-aggregate(phyl.mat~metadat$Fraction,FUN=mean)

rownames(phyl.ag)<-phyl.ag[,1]
phyl.ag.2<-as.matrix(phyl.ag[,-1])

phyl.ag.names<-rownames(phyl.ag)

windows()
barplot(apply(t(phyl.ag.2),2,rev),names=phyl.ag.names,col=c("#003f5c", "#444e86", "#955196","#dd5182","#ff6e54", "#ffa600"),horiz=FALSE,
        ylab="Relative abundance (%)")
# data wrangling  
phyl.df<- as.data.frame(phyl.ag.2)
phyl.df[,8]<-row.names(phyl.ag.2)
phy.df.t<-gather(phyl.df, "taxa",value, 1:7)


# Stacked + percent
windows()
ggplot(phy.df.t, aes(fill=taxa, y=value, x=V8)) + 
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  ggtitle("Top Phylum") +
  theme_bw(base_size = 14) +
  xlab("site")+
  ylab("relative abundance %")


# 


