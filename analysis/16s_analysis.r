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


################## Load required libraries ### ############

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
library(readxl)
#install.packages(ggtree)
#install.packages("TreeTools")
library('TreeTools')

#-----load functions-------------#

# load in the function for making a heatmap with the tree #
heatmap.phylo <- function(x, Rowp, Colp, ...) {
  l = length(seq(.1, 6, .1))
  pal = colorRampPalette(c("#fffaa2", "#bb0000"))(l)
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
        xaxs="i", yaxs="i", axes=FALSE, xlab= "", ylab= "", breaks=seq(.1,6.1,.1))
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", yaxs="i", xlim=c(0,2), ylim=yl)
  text(rep(0,nrow(x)),1:nrow(x), row_order, pos=4,
       cex=1, xpd=NA)
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", xaxs="i", ylim=c(0,2), xlim=xl)
  text(1:ncol(x),rep(2,ncol(x)), col_order, srt=90, adj=c(1,.5),
       cex=1.5)
}

### load colors 
# colors :)
gold <- "#FFB000"
lightgold <- "#fffaa2"
purple <- "#785EF0"
blue <- "#739AFF"
pink <- "#DC267F"
lightpink <- "#ffc6e5"
lightpurple <- "#8563C1"

# poster colors 
#navyblue<-"#161D64"
#dustypurple<-"#A79EB6"

#purple<-"#543A81"
#lightblue <-"#739AFF"
#blue<-"#5560A7"
#pink <- "#E48888"
#lightpink<- "#FEDDD2"
#gold<- "#AE330D"
#lightgold<-"#E48888"



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
#asvs.phyloseq<- (asvs.t)
#taxon<-taxon[,1:7]
#metadat<-as.matrix(metadat)
#y<-colnames(asvs.raw)
#rownames(metadat) <- y
#metadat<-as.data.frame(metadat)

#import it phyloseq
#Workshop_ASVS <- otu_table(asvs.phyloseq, taxa_are_rows = FALSE)
#Workshop_metadat <- sample_data(metadat)
#Workshop_taxo <- tax_table(as.matrix(taxon))
#ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

# take a look at PS object
#print(ps)
# we a get a Plyloseq object with  15027 taxa

## Determine minimum available reads per sample ##
min.s<-min(rowSums(asvs.t))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs.t, min.s)
#------------ rarefaction curve------------#
#S <- specnumber(asvs.t) # observed number of species
#raremax <- min(rowSums(asvs.t))
#plot(asvs.t, asvs.r, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
#abline(0, 1)
#out<-rarecurve(asvs.t, step = 20, sample = raremax, col = "blue", cex = 0.6)
## build plots
#col <- c("black")
#lty <- c("solid")
#lwd <- rep(1, 42)
#pars <- expand.grid(col = col, lty = lty, lwd = lwd, 
#                    stringsAsFactors = FALSE)
#Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
#Smax <- sapply(out, max)
## save
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/rarecurve.svg",width = 6, height=6 )
#plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Sequences",
#     ylab = "ASV", type = "n")
#abline(v = raremax)
#for (i in seq_along(out)) {
#  N <- attr(out[[i]], "Subsample")
#  with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
#}
#dev.off()


###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

##------make phyloseq object with rarefied data -------#

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
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

#test it worked
#sample_names(ps)
print(ps)
# 15027 taxa

# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that aren; in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
# 14593 taxa

#asvs.clean<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))



# rareferying by number of cells 
## Set the working directory; modify to your own ###
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###


#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/sequencing")
#n_cells <- read_excel("Cell_counts_For_sequencing.xlsx", na = "NA")
#n_cells<-select(n_cells, Sample_ID, Number_cells)
#n_cells

#data wrangleing and importing into Phyloseq#
## recode metadata 
#metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
#metadat<-metadat[, c(1,3:6)]
#metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
#metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

## Transpose ASVS table ##
#otu.t <- t(otu.raw)
## order metadata
#metadat<-metadat[order(metadat$SampleID),]
## order asvs table
#otu.t<-otu.t[order(row.names(otu.t)),]

## rearrange data for phyloseq
#otu.phyloseq<- (otu.t)
#taxon<-taxon[,1:7]
#metadat<-as.matrix(metadat)
#y<-colnames(otu.raw)
#rownames(metadat) <- y
#metadat<-as.data.frame(metadat)
## import it phyloseq
#Workshop_OTU <- otu_table(otu.phyloseq, taxa_are_rows = FALSE)
#Workshop_metadat <- sample_data(metadat)
#Workshop_taxo <- tax_table(as.matrix(taxon))
#ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
#ps
# 7727 taxa
# grab samples of total cells and boncat pos
#ps1 <- subset_samples(ps,Fraction =="BONCAT_Active"  |Fraction=="Total_Cells")
#ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
#any(taxa_sums(ps1) == 0)
#ps1
# 3208 taxa
#sample_data(ps1)

#normalize by the number of cells in each sample.
#otu.cells<-as.data.frame(otu_table(ps1))
#n_cells<-na.omit(n_cells)
#n_cells

#otu.cells<-cbind(otu.cells, n_cells)

#seq<-rowSums(otu.cells)
#otu.cells<-otu.cells*n_cells$Number_cells/seq

### looks great now put it back into phyloseq.

## rearrange data for phyloseq
#otu.phyloseq<- (otu.cells)

## import it phyloseq
#Workshop_OTU <- otu_table(otu.phyloseq, taxa_are_rows = FALSE)
#Workshop_metadat <- sample_data(metadat)
#Workshop_taxo <- tax_table(as.matrix(taxon))
#ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)
#ps
# 3208 taxa


#######------ remove chloroplasts-------#
# remove chloroplast DNA
#ps<-subset_taxa(ps, Class!=" Chloroplast")
#ps<-subset_taxa(ps, Genus!=" Mitochondria")
#ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
#ps<-prune_taxa(taxa_sums(ps) > 0, ps)
#any(taxa_sums(ps) == 0)
#ps

# 3173 taxa
# output df 
#asvs.clean<-as.data.frame(t(as.data.frame(otu_table(ps))))
#taxon<-as.data.frame(tax_table(ps))



########------DIVERSITY  figure---------########
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
rich %>%
arrange( -Observed)

rich<- rich %>%  filter(Plant!="NOPLANT", Fraction!="Inactive")
rich$compartment_BCAT <-factor(rich$compartment_BCAT, levels = c("Bulk_SoilTotal_DNA", "RhizosphereTotal_DNA", "RhizosphereTotal_Cells", "RhizosphereBONCAT_Active",
                                                                 "RootsTotal_Cells"   ,  "RootsBONCAT_Active" , "NoduleTotal_Cells" ,  "NoduleBONCAT_Active" ))          

# total + active number otus
svg(file="figures/16s/observed_divserity.svg",width = 6, height=7 )
windows(width = 6, height=7)
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=Observed, fill = Compartment, col= Compartment))+
  geom_boxplot() +
  scale_colour_manual(values = c("black", "black", gold, pink))+
  scale_fill_manual( values = c(blue, purple, lightgold, lightpink))+
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.position =  c(0.8, 0.8))+
  ylab("Number of ASVS")+
  xlab("Compartment")
dev.off()

## just active alpha
svg(file="figures/16s/observed_divserity_act.svg",width = 6, height=7 )
windows(width = 6, height=7)
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", Fraction=="BONCAT_Active")%>%
  ggplot(aes(x=compartment_BCAT, y=Observed, fill = Compartment, col= Compartment))+
  geom_boxplot() +
  scale_colour_manual(values = c( "black", gold, pink))+
  scale_fill_manual( values = c( purple, lightgold, lightpink))+
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.position =  c(0.8, 0.8))+
  ylab("Number of ASVS")+
  xlab("Compartment")
dev.off()




# total + active shannon
svg(file="figures/16s/shannon_diversity.svg",width = 6, height=7 )
windows(width = 6, height=7)
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=Shannon, fill = Compartment, col= Compartment))+
  geom_boxplot() +
  scale_colour_manual(values = c( blue,  purple, gold, pink))+
  scale_fill_manual( values = c(blue, lightpurple, lightgold, lightpink))+
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.position =  c(0.8, 0.8))+
  ylab("Shannon Diversity")+
  xlab("Compartment")
dev.off()



### t test

df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)

m1<-lm(Observed ~ Compartment+ BONCAT + Compartment*BONCAT,  data = df)
summary(m1)

#rhizosphere verse nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule")

m1<-lm(Observed ~ Compartment+ BONCAT + Compartment*BONCAT,  data = df)
summary(m1)

df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="POS")

m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)
# compartments are different 
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="SYBR")

m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)
# compartments are different 




#rhizosphere verse roots
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment== "Roots")

m1<-lm(Observed ~ Compartment+ BONCAT + Compartment*BONCAT,  data = df)
summary(m1)

df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots", BONCAT=="POS")

m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)
# compartments are different 
#CompartmentRoots  -895.35      45.35  -19.74 2.14e-07 ***

# rhizo sybr verse endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots", BONCAT=="SYBR")

m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)

  

# in nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule")
m1<-lm(Observed ~ BONCAT,  data = df)
summary(m1)
### boncat pos and active are different 
### BONCATSYBR    12.350      4.716   2.619  0.03447 * 

# in endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots")
m1<-lm(Observed ~ BONCAT,  data = df)
summary(m1)
### boncat pos and active are different 
### BONCATSYBR   -54.150     13.032  -4.155  0.00427 ** 

# in rhizosphere sybr verse boncat
df<-rich %>%
filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="DNA") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Observed ~ BONCAT,  data = df)
summary(m1)
# sybr and boncat pos are different
#BONCATSYBR    316.65      68.69    4.61  0.00246 ** 
  

# in rhizosphere DNA verse sYBR
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="POS") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Observed ~ BONCAT,  data = df)
summary(m1)
# sybr and DNA are different
#BONCATSYBR   -267.00      59.18  -4.512  0.00197 ** 
  

# in rhizosphere DNA verse boncat
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="SYBR") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Observed ~ BONCAT,  data = df)
summary(m1)
# sybr and DNA are different
#BONCATSYBR   -267.00      59.18  -4.512  0.00197 ** 


# bulk soil dna verse rhizosphere DNA
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="DNA") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Bulk_Soil")
m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)
## rhizosphere DNA verse bulk soil dna are not different
#CompartmentRhizosphere  -100.60      96.51  -1.042     0.32    



# nodule verse endosphere boncart
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="POS") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule"| Compartment=="Roots")
m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)

# nodule verse endosphere total

df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule"| Compartment=="Roots")
m1<-lm(Observed ~ Compartment,  data = df)
summary(m1)






#####stats for diversity #########
### t test
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)

m1<-lm(Shannon ~ Compartment+ BONCAT + Compartment*BONCAT,  data = df)
summary(m1)

#rhizosphere verse nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule")

m1<-lm(Shannon ~ Compartment+ BONCAT + Compartment*BONCAT,  data = df)
summary(m1)

df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="POS")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)
# compartments are different 
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="SYBR")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)
# compartments are different 

# roots boncat active verse total nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(compartment_BCAT=="NoduleTotal_Cells"| compartment_BCAT=="RootsBONCAT_Active")

m1<-lm(Shannon ~ compartment_BCAT,  data = df)
summary(m1)
# compartments are different 



#rhizosphere verse roots
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment== "Roots")

m1<-lm(Shannon ~ Compartment+ BONCAT + Compartment*BONCAT,  data = df)
summary(m1)

df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots", BONCAT=="POS")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)
# compartments are different 
#CompartmentRoots  -895.35      45.35  -19.74 2.14e-07 ***

# rhizo sybr verse endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots", BONCAT=="SYBR")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)

# rhizo boncat verse endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots", BONCAT=="POS")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)

# rhizo boncat verse nod boncat
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="POS")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)

# rhizo boncat verse nod sybr
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="SYBR")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)

# in nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)
### boncat pos and active are different 
### BONCATSYBR    12.350      4.716   2.619  0.03447 * 

# in endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)


# in rhizosphere sybr verse boncat
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="DNA") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)
#  0.13285    0.10457   1.271    0.244    

# in rhizosphere DNA verse sYBR
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="POS") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)
# sybr and DNA are different


# in rhizosphere DNA verse boncat
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="SYBR") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)


# bulk soil dna verse rhizosphere DNA
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="DNA") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Bulk_Soil")
m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)



# nodule verse endosphere bcat pos
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="POS") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule"| Compartment=="Roots")
m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)

# nod vs endo sybr
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule"| Compartment=="Roots")
m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)







############ taxa exploration  ######################

######### remove rare taxa, is the endosphyte diversity still higher in active??

root_total<-subset_samples(ps, Compartment=="Roots")
sample_names(root_total)
root_total<-prune_taxa(taxa_sums(root_total) > 0, root_total)
any(taxa_sums(root_total) == 0)
root_total
#603 taxa

## fitler for taxa that are rare
root_total<-prune_taxa(taxa_sums(root_total) > 15, root_total)
any(taxa_sums(root_total) == 0)
root_total
# 186 taxa

x<-otu_table(root_total)

#calc richness
rich<-estimate_richness(root_total, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

# make little metadata
metadat2<-metadat %>% filter(Compartment == "Roots")

# Data wrangling fo rdiversity of active microbes in each fraction
rich<-cbind(rich, metadat2)
rich<-as.data.frame(rich)
colnames(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))
rich %>%
  arrange( -Observed)

## t test
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)

## t test
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots")
m1<-lm(Observed ~ BONCAT,  data = df)
summary(m1)

########### more exploration
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

n<-sample_names(root_total)
root_total<-subset_samples(ps, Compartment=="Roots")
root_total<-prune_taxa(taxa_sums(root_total) > 0, root_total)
any(taxa_sums(root_total) == 0)
root_total
#406 taxa between both active and inactive
otu<-as.data.table(otu_table(root_total))
taxon<-(tax_table(root_total))

otu<-t(otu)
otu
total<-cbind(otu.clean,taxon)
#sum<-rowSums(total[,1:4])
#total<-mutate(total, sum= sum)
window("total")
hist(total$sum)

####### not including taxa that are really rare less than 5 reads 
####### are all those taxa in the root present in the bulk soil? ##########
otu<-as.data.frame(otu)
colnames(otu)
otu1<-filter(otu, V1 > 5 | V2 > 5 | V3 > 5 | V4 > 5 | V5 > 5 | V6 > 5 | V7 > 5 | V8 > 5 )
dim(otu1)
# 194 taxa
r_otus<- row.names(otu1)
###### okay now grad the total dna list
sample_data(ps)
total<-subset_samples(ps, Fraction=="Total_DNA")
total<-prune_taxa(taxa_sums(total) > 0, total)
any(taxa_sums(total) == 0)
total<-as.data.frame(t(otu_table(total)))

# 6041 taxa
t_otus<-row.names(total)
length(intersect(t_otus, r_otus)) # Apply setdiff function to see what's missing from the total DNA
length(setdiff(r_otus, t_otus))# 
# 62 taxa missing
mynames<-setdiff(r_otus, t_otus)

62/(132+62)


mynames

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



#########################VENN DIAGRAM are most taxa generalists or specialists?#################################
df<-subset_samples(ps, Fraction=="BONCAT_Active")
df<-prune_taxa(taxa_sums(df) > 0, df)
# subset for active
df<-prune_taxa(taxa_sums(df) > 5, df)
any(taxa_sums(df) == 0)

# don't include taxa that are super rare less then 5 reads
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))

dim(df)
df
n<-row.names.data.frame(df)
n[grepl("R" , n)]="rhizo"
n[grepl("E" , n)]="roots"
n[grepl("N" , n)]="nodule"

# summ by group
df<-rowsum(df, n)

df<-t(df)
df<-as.data.frame(df)
dim(df)
df
habitat<-rep(NA, 2621)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>1]="rhizo specailist"
habitat[df$nodule==0 & df$roots>1 & df$rhizo==0]="root specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
habitat
unique(habitat)

# rhizo   root    nod
# 1       0         0     -> rhizo specialist
# 0       1         0     -> root specialist
# 0       0         1     -> nod specialist
# 1       1         0     -> rhizo and root generalist
# 0       1         1     -> plant generalist
# 1       1         1     -> hyper generalist

df_habitat<-cbind(df, habitat)
df_habitat
############### calucating sums by hand
df_sums<-rowsum(df, habitat)
df_sums$rhizo
value <- df_sums$rhizo
value[c(2,3)]<-df_sums$nodule[c(2,3)]
value[c(6,7)]<-df_sums$roots[c(6,7)]
value
df_sums<-mutate(df_sums, value = value)

################ df wrangling for venndriagram
df[df>1] <- 1

nodule<-rownames(df[df$nodule==1,])
roots<-rownames(df[df$roots==1,])
rhizo<-rownames(df[df$rhizo==1,])

x <- list(
  nodule = nodule, 
  roots = roots, 
  rhizo = rhizo
)
x

################################## venn diagram #
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
ggvenn(
  x, 
  fill_color = c(pink, "#EFC000FF", purple),
  stroke_size = 0.5, set_name_size = 5,
  auto_scale = TRUE
)

####genus level??




###########-------- total cells 
df<-subset_samples(ps, Fraction=="Total_Cells")
df<-prune_taxa(taxa_sums(df) > 0, df)
# subset for active
#df<-prune_taxa(taxa_sums(df) > 5, df)
any(taxa_sums(df) == 0)

# don't include taxa that are super rare less then 5 reads
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))

dim(df)
head(df)
n<-row.names.data.frame(df)
n[grepl("N" , n)]="nodule"
n[grepl("E" , n)]="roots"
n[grepl("R" , n)]="rhizo"
n

# summ by group
df<-rowsum(df, n)
dim(df)
df<-t(df)
df<-as.data.frame(df)
dim(df)
head(df)
#habitat<-rep(NA, 2621)

#habitat[df$nodule==0 & df$roots==0 & df$rhizo>1]="rhizo specailist"
#habitat[df$nodule==0 & df$roots>1 & df$rhizo==0]="root specialist"
#habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
#habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
#habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
#habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
#habitat
#unique(habitat)

# rhizo   root    nod
# 1       0         0     -> rhizo specialist
# 0       1         0     -> root specialist
# 0       0         1     -> nod specialist
# 1       1         0     -> rhizo and root generalist
# 0       1         1     -> plant generalist
# 1       1         1     -> hyper generalist

#df_habitat<-cbind(df, habitat)
#df_habitat
############### calucating sums by hand
#df_sums<-rowsum(df, habitat)
#df_sums$rhizo
#value <- df_sums$rhizo
#value[c(2,3)]<-df_sums$nodule[c(2,3)]
#value[c(6,7)]<-df_sums$roots[c(6,7)]
#value
#df_sums<-mutate(df_sums, value = value)

################ df wrangling for venndriagram
df[df>1] <- 1

nodule<-rownames(df[df$nodule==1,])
roots<-rownames(df[df$roots==1,])
rhizo<-rownames(df[df$rhizo==1,])

x <- list(
  nodule = nodule, 
  roots = roots, 
  rhizo = rhizo
)
x

################################## venn diagram #
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
ggvenn(
  x, 
  fill_color = c(pink, "#EFC000FF", purple),
  stroke_size = 0.5, set_name_size = 5
)


####Bar chart - are taxa that are generalist realy abundant?

df<-subset_samples(ps, Fraction=="BONCAT_Active")
df<-prune_taxa(taxa_sums(df) > 0, df)
# subset for active
df<-prune_taxa(taxa_sums(df) > 5, df)
any(taxa_sums(df) == 0)

# don't include taxa that are super rare less then 5 reads
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))

dim(df)
df
n<-row.names.data.frame(df)
n[grepl("R" , n)]="rhizo"
n[grepl("E" , n)]="roots"
n[grepl("N" , n)]="nodule"

# summ by group
df<-rowsum(df, n)

df<-t(df)
df<-as.data.frame(df)
dim(df)
df
habitat<-rep(NA, 2621)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>1]="rhizo specailist"
habitat[df$nodule==0 & df$roots>1 & df$rhizo==0]="root specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="nodule and rhizo generalist"


habitat
unique(habitat)

# rhizo   root    nod
# 1       0         0     -> rhizo specialist
# 0       1         0     -> root specialist
# 0       0         1     -> nod specialist
# 1       1         0     -> rhizo and root generalist
# 0       1         1     -> plant generalist
# 1       1         1     -> hyper generalist



df_habitat<-cbind(df, habitat)
dim(df_habitat)
unique(df_habitat$habitat)

df_habitat$abundace<- rowSums(df_habitat %>% 
            dplyr::select(contains("O"))) %>% 
  glimpse()

df_habitat$mean_abundance<- rowMeans(df_habitat %>% 
             dplyr::select(contains("O"))) %>% 
  glimpse()

head(df_habitat)

df_habitat%>% filter(rhizo>1)%>%
  ggplot()+
  geom_boxplot(mapping=aes(x = reorder(habitat, +rhizo), y = log10(rhizo)))+
  #geom_jitter(mapping=aes(x = reorder(habitat, +rhizo), y = log10(rhizo)))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=0.9))+
  ylab("log 10 relative abundance in rhizosphere")

df_habitat%>% filter(roots>1)%>%
  ggplot()+
  geom_boxplot(mapping=aes(x = reorder(habitat, +roots), y = log10(roots)))+
  geom_jitter(mapping=aes(x = reorder(habitat, +roots), y = log10(roots)))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=0.9))+
  ylab("log 10 relative abundance in root")


df_habitat%>% filter(nodule>1)%>%
  ggplot()+
  geom_boxplot(mapping=aes(x = reorder(habitat, +nodule), y = log10(nodule)))+
  geom_jitter(mapping=aes(x = reorder(habitat, +nodule), y = log10(nodule)))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=0.9))+
  ylab("log 10 relative abundance in root")








############################PCOA plots  ########################
#Pcoa on rarefied asvs Data

# rm ctl
ps<-subset_samples(ps, Compartment !="ctl")

# Calculate Bray-Curtis distance between samples

otus.bray<-vegdist(otu_table(ps), method = "bray")
ps


# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(40-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
# swtich order
switch<-otus.p[,2:1]

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
metadat<-filter(metadat, Compartment!="ctl")
as.factor(metadat$Compartment)
as.factor(metadat$Fraction)
as.factor(metadat$compartment_BCAT)
levels(as.factor(metadat$compartment_BCAT))
levels(as.factor(metadat$Fraction))

# colors :)
#gold <- "#FFB000"
#purple <- "#785EF0"
#blue <- "#739AFF"
#pink <- "#DC267F"

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/pcoa")

svg(file="Pcoa_plantvsoil_raw.svg",width = 6, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 6, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, col=c("black"),
       pch=c(circle, triangle, diamond)[as.factor(metadat$Fraction)],
       lwd=1,cex=2,
       bg=c( blue, pink , pink, purple, purple, purple, gold, gold)[as.factor(metadat$compartment_BCAT)])
ordiellipse(otus.pcoa, metadat$Compartment,  
            kind = "se", conf=0.95, #label=T, 
            #draw = "polygon",
            lwd=2, col="black")


legend("top",legend=c("beads", "Bulk_SoilTotal_DNA" ,      "neg control"   ,  "NoduleBONCAT_Active"    ,  "Nodule Total Cells"     ,      "RhizosphereBONCAT_Active",
                      "Rhizosphere Total Cells","RhizosphereTotal_DNA" ,    "Roots BONCAT Active"   ,    "Roots Total Cells"),
       pch=c(15,5,0,1,2 , 1, 2 , 5, 1, 2),
       col= c("black", "#739AFF", "black", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000"),
       bty = "n",
       inset = c(.05, 0))

dev.off()



##########-------run a new pcoa on just soil <3
ps
metadat
ps2<-subset_samples(ps, Compartment !=  "Nodule" & Compartment != "Roots" & Compartment !="ctl" & Fraction != "Total_DNA")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 13390 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")


# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(8-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

# subset metadata1
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots" & Compartment!="ctl" & Fraction != "Total_DNA")

as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

square <- 22
diamond <- 23
triangle <- 24
circle <- 21

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/soil_raw.svg",width = 4, height=4 )
windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(circle,triangle)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c(purple)[as.factor(metadat2$Compartment)])

ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "se", conf=0.95, label=T, 
            #draw = "polygon",
            lwd=2, col="black")

legend("center",legend=c( "Rhizosphere BONCAT_Active", "Rhizosphere Total Cells"), 
       pch=c(circle, triangle),
       cex=1.1, 
       fill=c(blue,"#739AFF"),
       bty = "n")

legend("topright", inset=c(-0.2,0), legend=c("Rhizosphere BONCAT_Active", "Rhizosphere Total Cells"), pch=c(circle, triangle))



dev.off()

####### active verse total RDA





#################-------------------just plant##
ps
ps2<-subset_samples(ps, Compartment !=  "Bulk_Soil" & Compartment != "Rhizosphere" & Compartment != "ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
df<-as.data.frame(otu_table(ps2))
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

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/plant_raw.svg",width = 4, height=4 )
windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(circle, triangle)[as.factor(metadat2$Fraction)],
       lwd=1,cex=2,
       bg=c(pink, pink, gold, gold)[as.factor(metadat2$compartment_BCAT)])
ordiellipse(otus.pcoa, metadat2$compartment_BCAT,  
            kind = "se", conf=0.95, label=T, 
            draw = "polygon",
            lwd=2, col=NULL)




#legend("topleft",legend=c("Flow cyto control", "Bulk soil total DNA", " PCR control",  "Rhizosphere BONCAT_Active" , "Rhizosphere Inactive",  "Rhizosphere Total DNA"), 
#       pch=c(15,5, 0, 1,2,5),
#       cex=1.1, 
#       col=c("black", "#739AFF",  "black", "#785EF0", "#785EF0",  "#785EF0"),
#       bty = "n")

dev.off()

#########-------------- permanova----------------#########

#full model
# use output that doesn't have chloroplast for permanova. 
dim(asvs.clean)
asvs.clean
otu.perm<- adonis2(t(asvs.clean)~ Compartment+ Fraction+Compartment*Fraction, data = metadat, permutations = 999, method="bray")
otu.perm
#Df SumOfSqs      R2       F Pr(>F)    
#Compartment           3   7.9377 0.68157 33.0965  0.001 ***
#  Fraction              2   0.8828 0.07580  5.5211  0.001 ***
#  Compartment:Fraction  2   0.2675 0.02297  1.6728  0.122    

#analysis of similarities
otu.ano<- anosim(t(otu.r.clean), grouping =  metadat$Compartment, permutations = 999)
summary(otu.ano)

#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=metadat$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse



##### subset data by compartment
# with in a compartment (rhizo, root, nodule) are the fractions evenly dispersed?
# just soil
ps
metadat
ps2<-subset_samples(ps, Compartment !=  "Nodule" & Compartment != "Roots" & Compartment !="ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 13390 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots" & Compartment!="ctl") 
#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$compartment_BCAT))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# unequal varience

# just rhizosphere
ps
metadat
ps2<-subset_samples(ps, Compartment ==  "Rhizosphere" & Compartment !="ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 8665 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment ==  "Rhizosphere" & Compartment!="ctl") 
#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# unequal variance


# just rhizosphere total cells and active
ps
metadat
ps2<-subset_samples(ps, Compartment ==  "Rhizosphere" & Compartment !="ctl" & Fraction!="Total_DNA")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 8665 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment ==  "Rhizosphere" & Compartment!="ctl" & Fraction!="Total_DNA")
#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# unequal variance


# just Total DNA
ps
metadat
ps2<-subset_samples(ps, Fraction=="Total_DNA" & Plant =="Clover")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 8665 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Fraction=="Total_DNA" & Plant=="Clover")
#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Compartment))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# unequal variance











######## Permanova just plant #############
ps
metadat
ps2<-subset_samples(ps, Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 691 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# equal variance :)
otu.perm<- adonis2((otu_table(ps2))~ Compartment+ Fraction + Compartment*Fraction, data = metadat2, permutations = 999, method="bray")
otu.perm
# significant effect of compartment and fraction
#                     Df SumOfSqs      R2       F Pr(>F)    
#Compartment           1  0.18467 0.36664 17.4150  0.001 ***
#  Fraction              1  0.13593 0.26986 12.8180  0.002 ** 
#  Compartment:Fraction  1  0.03463 0.06876  3.2659  0.076 .  
#Residual             14  0.14846 0.29474                   
#Total                17  0.50369 1.00000

### nodule active verse bulk
ps2<-subset_samples(ps, Compartment ==  "Nodule" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
# 131 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment ==  "Nodule" )
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# equal variance :)
otu.perm<- adonis2((otu_table(ps2))~ Fraction, data = metadat2, permutations = 999, method="bray")
otu.perm
# significant effect of compartment and fraction
#Fraction  1 0.035435 0.43128 5.3083  0.032 *
  

### endo active verse bulk
ps2<-subset_samples(ps, Compartment ==  "Roots" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
# 603 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment ==  "Roots" )
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
# equal variance :)
otu.perm<- adonis2((otu_table(ps2))~ Fraction, data = metadat2, permutations = 999, method="bray")
otu.perm
# significant effect of compartment and fraction
#Fraction  1  0.13512 0.57049 9.2975  0.012 *


# total cells endo vs nod

ps2<-subset_samples(ps, Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl" &Fraction ==  "Total_Cells" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 266 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl" &Fraction ==  "Total_Cells")


otu.perm<- adonis2((otu_table(ps2))~ Compartment, data = metadat2, permutations = 999, method="bray")
otu.perm
# significan

# boncat endo vs nod

ps2<-subset_samples(ps, Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl" &Fraction ==  "BONCAT_Active" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
# 266 asvs
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl" &Fraction == "BONCAT_Active")


otu.perm<- adonis2((otu_table(ps2))~ Compartment, data = metadat2, permutations = 999, method="bray")
otu.perm
#0.179721 0.67208 16.396  0.008 **


##rhizo active verse rhizo total cells
ps.r
ps2<-subset_samples(ps.r, Compartment ==  "Rhizosphere"  & Fraction!="Total_DNA")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
df<-as.data.frame(otu_table(ps2))
row.names(df)


metadat_y<-metadat[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_DNA") ,]
otu.perm<- adonis2(df~ Fraction, data = metadat_y, permutations = 999, method="bray")
otu.perm


dispersion <- betadisper(otus.bray, group=metadat$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse




#rhizo active ver total DNA
test<-otu.r.clean.t[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive"),]
metadat_y<-metadat[which(metadat$Compartment == "Rhizosphere" & metadat$Fraction != "Total_Cells" & metadat$Fraction != "Inactive") ,]
otu.perm<- adonis2(test~Fraction, data = metadat_y, permutations = 999, method = "bray")
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

#####-------PHYLOGENIC SIGNAL STARTS HERE ---------##########

## select active taxa
ps1<- subset_samples(ps, Fraction=="BONCAT_Active" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 10, ps1)
any(taxa_sums(ps1) == 0)
ps1
taxon<-as.data.frame(tax_table(ps1))

# 2169 taxa from rarefied ASVS
otu<-as.data.frame(t(as.data.frame(otu_table(ps1))))
dim(otu)
head(otu)
####### adding taxa information #######
# add taxa names
n<-row.names(otu)
otu<-mutate(otu, OTU=n)
n<-row.names(taxon)
taxon<- mutate(taxon, OTU=n)
otu<-left_join(otu, taxon)

# put g__ in front
otu$Genus <- gsub('\\ ' , '', otu$Genus)
otu$Genus<-paste0( "g__",otu$Genus)
otu$Family <- gsub('\\ ' , '', otu$Family)
otu$Family<-paste0( "f__",otu$Family)
otu$Order <- gsub('\\ ' , '', otu$Order)
otu$Order<-paste0( "o__",otu$Order)
otu$Class <- gsub('\\ ' , '', otu$Class)
otu$Class<-paste0( "c__",otu$Class)
otu$Phyla <- gsub('\\ ' , '', otu$Phyla)
otu$Phyla<-paste0( "p__",otu$Phyla)


# if genus is missing put family name, if still missing put order, class, domain, etc. 
otu$Genus[which(otu$Genus == 'g__')] <- otu$Family[which(otu$Genus == 'g__')]
otu$Genus[which(otu$Genus == 'f__')] <- otu$Order[which(otu$Genus == 'f__')]
otu$Genus[which(otu$Genus == 'o__')] <- otu$Class[which(otu$Genus == 'o__')]
otu$Genus[which(otu$Genus == 'c__')] <- otu$Phyla[which(otu$Genus == 'c__')]
otu$Genus[which(otu$Genus == 'p__')] <- otu$Domain[which(otu$Genus == 'p__')]

# separate archaea and bacteria
#otu<-otu %>% filter( Domain == "Bacteria")
#otu_arch<-otu %>% filter( Domain == "Archaea")


# aggregate by genus level 
#otu<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Genus, data = otu, FUN = sum, na.rm = TRUE)
#otu_arch<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Genus, data = otu_arch, FUN = sum, na.rm = TRUE)

# make row names genuses
#row.names(otu) <- otu$Genus
#row.names(otu_arch) <- otu_arch$Genus

# how many genuses are there?
length(unique(otu$Genus))
# 395 bacteria unique taxa groups


###################changing my data's names so they match the tree ####################

# put order as same for type III 
otu$Genus[grepl("g__type_III", otu$Genus )] <- otu$Order[grepl("g__type_III", otu$Genus) ]


#g__Tepidisphaeraceae" is a family
otu$Genus[grepl("g__Tepidisphaeraceae", otu$Genus)] <- otu$Family[grepl("g__Tepidisphaeraceae", otu$Genus)]

# switch "g__Nitrospira_ to "g__Nitrospira_A" 
otu$Genus <- sub("g__Nitrospira" ,  "g__Nitrospira_A"   ,  otu$Genus )

# "g__Vicinamibacteraceae" is a family
otu$Genus[grepl("g__Vicinamibacteraceae", otu$Genus)] <-  otu$Family[grepl("g__Vicinamibacteraceae", otu$Genus )] 

#  "g__Vampirovibrionaceae"  is a family
otu$Genus[grepl("g__Vampirovibrionaceae" , otu$Genus)] <-  otu$Family[grepl("g__Vampirovibrionaceae" , otu$Genus )] 

# "g__Vermiphilaceae"     is family
otu$Genus[grepl( "g__Vermiphilaceae"   , otu$Genus)] <-  otu$Family[grepl(  "g__Vermiphilaceae" , otu$Genus )] 

# "g__Vampirovibrionales"  is a order
otu$Genus[grepl( "g__Vampirovibrionales" , otu$Genus)] <-  otu$Order[grepl("g__Vampirovibrionales" , otu$Genus )] 


#### change g__Peribacteria" to  "g__Peribacter"
otu$Genus <- sub("g__Peribacteria", "g__Peribacter" ,  otu$Genus )

# Parcubacteria is a phyla so make it a phyla
otu$Genus[grepl("g__Parcubacteria", otu$Genus)] <-  otu$Phyla[grepl("g__Parcubacteria", otu$Genus )] 


### _Pedosphaeraceae is a family, call it a family
otu$Genus[grepl("g__Pedosphaeraceae", otu$Genus)] <-  otu$Family[grepl("g__Pedosphaeraceae", otu$Genus )] 


### replace unknown family with order
otu$Genus[grepl("g__Unknown_Family", otu$Genus)] <-  otu$Order[grepl("g__Unknown_Family", otu$Genus )] 

### Berkelbacteria is actual a phyla
otu$Genus[grepl("g__Berkelbacteria", otu$Genus)] <-  otu$Phyla[grepl("g__Berkelbacteria", otu$Genus )] 

### Gracilibacteria is actual a phyla
otu$Genus[grepl("g__Gracilibacteria", otu$Genus)] <-  otu$Phyla[grepl("g__Gracilibacteria", otu$Genus )] 

#"g__Chthoniobacteraceae" is an family, but only 1 genus in it, so I'm gonna call it that 
otu$Genus <- sub("g__Chthoniobacteraceae", "g__Chthonomonas" ,  otu$Genus , ignore.case = TRUE)

#"g__Chthoniobacteraceae" is an order, but only 1 genus in it, so I'm gonna call it that 
otu$Genus <- sub("g__Chthonomonadales", "g__Chthonomonas" ,  otu$Genus , ignore.case = TRUE)

#"g__#Entotheonaella" is an family, but only 1 genus in it, so I'm gonna call it that 
otu$Genus <- sub("g__Entotheonellaceae", "g__Entotheonella" ,  otu$Genus )

#Fimbriimonadaceae is a family  , but only 1 genus in it, so I'm gonna call it that 
otu$Genus <- sub("g__Fimbriimonadaceae", "g__Fimbriimonas" ,  otu$Genus )

# Fimbriimonadales is an order , but only 1 genus in it, so I'm gonna call it that 
otu$Genus <- sub("g__Fimbriimonadales", "g__Fimbriimonas" ,  otu$Genus )

# Hyphomicrobium call it g__Hyphomicrobium_A
otu$Genus <- sub("g__Hyphomicrobium", "g__Hyphomicrobium_A" ,  otu$Genus )

### change "g__Latescibacteraceae" to lascibacter
otu$Genus <- sub("g__Latescibacterota", "g__Latescibacter" ,  otu$Genus )

### change "g__Latescibacterota" to lascibacter
otu$Genus <- sub("g__Latescibacteraceae", "g__Latescibacter" ,  otu$Genus )

## change  to g__Massilia to "g__Massilibacterium"  because they are the same
otu$Genus <- sub("g__Massilia", "g__Massilibacterium" ,  otu$Genus )

# Methyloligellaceae is a family so call it a family
otu$Genus[grepl("Methyloligellaceae", otu$Genus)] <-  otu$Family[grep("Methyloligellaceae", otu$Genus )] 

# Microgenomatia is a phyla so call it a phyla
otu$Genus[grepl("g__Microgenomatia", otu$Genus)] <-  otu$Phyla[grep("g__Microgenomatia", otu$Genus )] 

# Rokubacteriales is an order. call it an order
otu$Genus[grepl("g__Rokubacteriales", otu$Genus)] <-  otu$Order[grep("g__Rokubacteriales", otu$Genus )] 


# replace linage with phyla
otu$Genus[grepl("g__Lineage", otu$Genus)] <-  otu$Phyla[grep("g__Lineage", otu$Genus )] 

## anything with number replace with a higher taxonomic level
otu$Genus[grep("[7]", otu$Genus)] <-  otu$Family[grep("[7]", otu$Genus )] 
otu$Genus[grep("[7]", otu$Genus)] <-  otu$Order[grep("[7]", otu$Genus )] 
otu$Genus[grep("[7]", otu$Genus)] <-  otu$Class[grep("[7]", otu$Genus )] 
otu$Genus[grep("[7]", otu$Genus)] <-  otu$Phyla[grep("[7]", otu$Genus )] 

otu$Genus[grep("[-]", otu$Genus)] <-  otu$Family[grep("[-]", otu$Genus )] 
otu$Genus[grep("[-]", otu$Genus)] <-  otu$Order[grep("[-]", otu$Genus )] 
otu$Genus[grep("[-]", otu$Genus)] <-  otu$Class[grep("[-]", otu$Genus )] 
otu$Genus[grep("[-]", otu$Genus)] <-  otu$Phyla[grep("[-]", otu$Genus )] 

otu$Genus[grep("[1]", otu$Genus)] <-  otu$Family[grep("[1]", otu$Genus )] 
otu$Genus[grep("[1]", otu$Genus)] <-  otu$Order[grep("[1]", otu$Genus )] # A0
otu$Genus[grep("[1]", otu$Genus)] <-  otu$Class[grep("[1]", otu$Genus )] 
otu$Genus[grep("[1]", otu$Genus)] <-  otu$Phyla[grep("[1]", otu$Genus )] 
otu$Genus[grep("[1]", otu$Genus)] <-  otu$Domain[grep("[1]", otu$Genus )] 

otu$Genus[grep("[2]", otu$Genus)] <-  otu$Family[grep("[2]", otu$Genus )] 
otu$Genus[grep("[2]", otu$Genus)] <-  otu$Order[grep("[2]", otu$Genus )] 
otu$Genus[grep("[2]", otu$Genus)] <-  otu$Class[grep("[2]", otu$Genus )] 
otu$Genus[grep("[2]", otu$Genus)] <-  otu$Phyla[grep("[2]", otu$Genus )] 
otu$Genus[grep("[2]", otu$Genus)] <-  otu$Domain[grep("[2]", otu$Genus )] 

otu$Genus[grep("[3]", otu$Genus)] <-  otu$Family[grep("[3]", otu$Genus )] 
otu$Genus[grep("[3]", otu$Genus)] <-  otu$Order[grep("[3]", otu$Genus )] 

otu$Genus[grep("[4]", otu$Genus)] <-  otu$Family[grep("[4]", otu$Genus )] 
otu$Genus[grep("[4]", otu$Genus)] <-  otu$Order[grep("[4]", otu$Genus )] 
otu$Genus[grep("[4]", otu$Genus)] <-  otu$Class[grep("[4]", otu$Genus )] 
otu$Genus[grep("[4]", otu$Genus)] <-  otu$Phyla[grep("[4]", otu$Genus )] 
otu$Genus[grep("[4]", otu$Genus)] <-  otu$Domain[grep("[4]", otu$Genus )] 

otu$Genus[grep("[5]", otu$Genus)] <-  otu$Family[grep("[5]", otu$Genus )] 
otu$Genus[grep("[5]", otu$Genus)] <-  otu$Order[grep("[5]", otu$Genus )] 
otu$Genus[grep("[5]", otu$Genus)] <-  otu$Class[grep("[5]", otu$Genus )] 
otu$Genus[grep("[5]", otu$Genus)] <-  otu$Phyla[grep("[5]", otu$Genus )] 
otu$Genus[grep("[5]", otu$Genus)] <-  otu$Domain[grep("[5]", otu$Genus )] 

otu$Genus[grep("[6]", otu$Genus)] <-  otu$Family[grep("[6]", otu$Genus )] 

otu$Genus[grep("[8]", otu$Genus)] <-  otu$Order[grep("[8]", otu$Genus )] 
otu$Genus[grep("[8]", otu$Genus)] <-  otu$Class[grep("[8]", otu$Genus )]

otu$Genus[grep("[9]", otu$Genus)] <-  otu$Order[grep("[9]", otu$Genus )] 


# swap "o__Absconditabacteriales_(SR1)" for "g__Absconditicoccus"]
otu$Genus <- sub("\\o__Absconditabacteriales_\\(SR1\\)", "o__Absconditabacteriales",  otu$Genus )

# swap "g__Clostridium_sensu_stricto_13"for "g__Clostridium" 
otu$Genus <- sub("g__Clostridium_sensu_stricto_13", "g__Clostridium" ,  otu$Genus , ignore.case = TRUE)


# candidatus, because that just means it hasn't been cultured
otu$Genus<- gsub('g__Candidatus_', 'g__', otu$Genus )

#change Armatimonadota to  "g__Armatimonas
otu$Genus<- gsub('g__Armatimonadales', 'g__Armatimonas', otu$Genus )

## change g__Altererythrobacter to g__Altererythrobacter_D" 
otu$Genus <- gsub('g__Altererythrobacter', 'g__Altererythrobacter_D', otu$Genus )

# for the groups with genusus that are just number use higher taxonomic classification
otu$Genus[grepl('g__0', otu$Genus )] <- otu$Class[grepl('g__0', otu$Genus)]

# for the groups with uncultured as genus  use higher taxonomic classification
otu$Genus[grepl('g__un', otu$Genus )] <- otu$Family[grepl('g__un', otu$Genus)]
otu$Genus[grepl('f__un', otu$Genus )] <- otu$Order[grepl('f__un', otu$Genus)]
otu$Genus[grepl('o__un', otu$Genus )] <- otu$Class[grepl('o__un', otu$Genus)]

# for the groups with subgroup as genus  use higher taxonomic classification
otu$Genus[grepl('Subgroup', otu$Genus )] <- otu$Family[grepl('Subgroup', otu$Genus)]
otu$Genus[grepl('f__Sub', otu$Genus )] <- otu$Order[grepl('f__Sub', otu$Genus)]
otu$Genus[grepl('o__Sub', otu$Genus )] <- otu$Class[grepl('o__Sub', otu$Genus)]

# rename rhizobium group as such
otu$Genus <- gsub('g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' , 'g__Rhizobium', otu$Genus)

# rename burkholderia group as such
otu$Genus <- gsub( 'g__Burkholderia-Caballeronia-Paraburkholderia' , 'g__Burkholderia', otu$Genus)

# replace Acidibacter with the order
otu$Genus <- gsub('g__Acidibacter',  'o__Gammaproteobacteria_Incertae_Sedis' , otu$Genus)

# replace "g__Absconditabacteriales_(SR1)") with the order level
otu$Genus[grepl('g__Abscond', otu$Genus )] <- otu$Order[grepl('Abscond', otu$Genus)]

################### aggregate by genus  #############

# how many genuses are there now?
length(unique(otu$Genus))
# 352 bacteria unique taxa groups

# aggregate by genus level 
otu<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Genus, data = otu, FUN = sum, na.rm = TRUE)

# make row names genuses
row.names(otu) <- otu$Genus
dim(otu)

# 352 genuses

# rename#  columns
#otu<-select(otu, -Phyla, -Domain, -Class, -Order, -Family, -Species, -OTU)
#n<-c("Genus", "Root4", "Nodule4", "Root1",  "Nodule1",  "Root2",  "Nodule2", "Root3", "Nodule3", "Rhizosphere4", "Rhizosphere1", "Rhizosphere2", "Rhizosphere3") 
#colnames(otu)<-n

## filter for taxa that were not present in the plant.
#otu1<-filter(otu, Root4 > 0 | Nodule4 > 0 | Root1 > 0 | Nodule1 > 0 | Root2 > 0 | Nodule2 > 0 | Root3 > 0 | Nodule3 > 0  )

#------> skip this part if you jsut want to group without phylogeny
 
#######################Read in the tree file and taxonomy file 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
#bacteria
#tree = read.tree("bac120.tree")
tree = read.tree("genus.nwk")
#tax = read.table("bac120_taxonomy.tsv")
map = read.table("family.map")




################# making tip names genus  ############ 
# put a g__ in front of all of these
tree$tip.label<-paste0("g__", tree$tip.label)
tree$tip.label
# remove random little /'-
tree$tip.label<-gsub("\\'", "", tree$tip.label)
length(tree$tip.label)
## 3474 genus

length(intersect(unique(otu$Genus), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
# 153 intersect
length(setdiff(unique(otu$Genus), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
# 199 don't match
mynames<-(setdiff(unique(otu$Genus), tree$tip.label))
mynames
#df<-otu %>% filter( Genus %in% mynames)

############ edit label in my data set so they match the df
###### exploration
mynames[grepl("g__", mynames)]
length(mynames[grepl("g__", mynames)])
# 59 genus
mynames[grepl("f__", mynames)]
length(mynames[grepl("f__", mynames)])
#76 families
mynames[grepl("o__", mynames)]
length(mynames[grepl("o__", mynames)])
#33 orders
mynames[grepl("c__", mynames)]
length(mynames[grepl("c__", mynames)])
# 20 classes 
mynames[grepl("p__", mynames)]
length(mynames[grepl("p__", mynames)])
# 11 phyla


### make a smaller tree for my data



#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, otu$Genus) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove,  collapse.singles = TRUE,) # remove asvs we don't need
length(tree.short$tip.label)
####### node labels
# add names of families
# grab names from map file
tree.label<-tree.short
(tree.label$node.label)
head(map$V1)
# remove the G0 ones
m<-map$V1[grepl("N", map$V1)]
map.short<-filter(map, map$V1 %in% m)
head(map.short)

# find nodes that match my data frame
x<-tree.label$node.label
df<-data.frame(x)
k<-intersect(df$x, map.short$V1)

#grab names
df1<-filter(map.short, V1 %in% k)
df1$V2
k
#df<-left_join(df, df1)
for (i in seq_along(k)) {
tree.label$node.label<-sub(k[i], df1$V2[i], tree.label$node.label)
}
tree.label$node.label

#### add names of orders
map = read.table("order.map")
# remove the G0 ones
#map
m<-map$V1[grepl("N", map$V1)]
map.short<-filter(map, map$V1 %in% m)
#dim(map.short)
# find nodes that match my data frame
x<-tree.label$node.label
df<-data.frame(x)
k<-intersect(df$x, map.short$V1)
k
#grab names
df1<-filter(map.short, V1 %in% k)
df1$V2
#k
#df<-left_join(df, df1)
for (i in seq_along(k)) {
  tree.label$node.label<-sub(k[i], df1$V2[i], tree.label$node.label)
}
tree.label$node.label

#### add names of class
map = read.table("class.map")
# remove the G0 ones
#map
m<-map$V1[grepl("N", map$V1)]
map.short<-filter(map, map$V1 %in% m)
#dim(map.short)

# find nodes that match my data frame
x<-tree.label$node.label
df<-data.frame(x)
k<-intersect(map.short$V1,df$x)
k
#grab names
df1<-filter(map.short, V1 %in% k)
df1$V2
#k
#df<-left_join(df, df1)
for (i in seq_along(k)) {
  tree.label$node.label<-sub(k[i], df1$V2[i], tree.label$node.label)
  
}
tree.label$node.label

#### add names of phyla
map = read.table("phylum.map")
# remove the G0 ones
#map
m<-map$V1[grepl("N", map$V1)]
map.short<-filter(map, map$V1 %in% m)
#dim(map.short)

# find nodes that match my data frame
x<-tree.label$node.label
df<-data.frame(x)
k<-intersect(map.short$V1,df$x)
k
#grab names
df1<-filter(map.short, V1 %in% k)
df1$V2
#k
#df<-left_join(df, df1)
for (i in seq_along(k)) {
  tree.label$node.label<-sub(k[i], df1$V2[i], tree.label$node.label)
  
}
tree.label$node.label


###### this didn't add all the names, but added some
tree.label$node.label<-sub("N14", "Bacteria", tree.label$node.label)
tree.label$node.label <- sub ("N1800", 	"Solirubrobacterales", tree.label$node.label)
tree.label$node.label <- gsub ("N3651", "Geodermatophilaceae", tree.label$node.label)
tree.label$node.label <- gsub ("N924", "Actinomycetota", tree.label$node.label)
tree.label$node.label <- gsub( "N3653", "Pseudonocardiaceae",  tree.label$node.label)
tree.label$node.label <- gsub( "N1301", "Armatimonadota",  tree.label$node.label)
tree.label$node.label <- gsub("N321", "Bacillota", tree.label$node.label)
tree.label$node.label <- gsub("N39", "Gram postive", tree.label$node.label)
tree.label$node.label <- gsub("N234", "Gram negative", tree.label$node.label)
tree.label$node.label <- gsub("N4739", "Actinomycetaceae", tree.label$node.label)
tree.label$node.label <- gsub("N2980", "Actinomycetia" , tree.label$node.label)
tree.label$node.label <- gsub("N4368", "Micrococcaceae", tree.label$node.label)
tree.label$node.label <- gsub("N1390",  "Pseudomonadota" , tree.label$node.label)
tree.label$node.label <- gsub("N2203" , "Betaproteobacteria", tree.label$node.label)
tree.label$node.label <- gsub("N3518", "Burkholderiales" , tree.label$node.label)
tree.label$node.label <- gsub("N8500", "Comamonadaceae", tree.label$node.label)
tree.label$node.label <- gsub("N3889", "Coxiellaceae", tree.label$node.label)
tree.label$node.label <- gsub("N2521", "Gammaproteobacteria", tree.label$node.label)
tree.label$node.label <- gsub("N5794", "Xanthomonadaceae", tree.label$node.label)
tree.label$node.label <- gsub("N4999", "superorder1", tree.label$node.label)
tree.label$node.label <- gsub("N3180", "superorder2", tree.label$node.label)
tree.label$node.label <- gsub("N7080", "subfamily1", tree.label$node.label)
tree.label$node.label <- gsub("N6230", "subfamily2", tree.label$node.label)


tree.label$node.label[grepl( "N31", tree.label$node.label)]

windows(10,10)
plot(tree.label, no.margin=TRUE,  cex = .5, show.node.label = TRUE)
nodelabels()
dev.off()
length(tree.label$tip.label)

tree.label$node.label


############## tips to add################

## add families
#mynames[grepl("f__", mynames)]
#            "f__Anaerolineaceae"            "f__Ardenticatenaceae"         
#[4] "f__Azospirillaceae"            "f__Bacillaceae"                "f__Bdellovibrionaceae"        
#[7] "f__Beijerinckiaceae"           "f__Blastocatellaceae"          "f__Burkholderiaceae"          
#[10] "f__Caldilineaceae"             "f__Caulobacteraceae"           "f__Chitinophagaceae"          
#[13] "f__Chloroflexaceae"            "f__Chthoniobacteraceae"        "f__Clostridiaceae"            
#[16] "f__Comamonadaceae"             "f__Devosiaceae"                "f__Diplorickettsiaceae"       
#[19] "f__Enterobacteriaceae"         "f__Erwiniaceae"                "f__Gemmataceae"               
#[22] "f__Gemmatimonadaceae"          "f__Geobacteraceae"             "f__Ilumatobacteraceae"        
#[25] "f__Intrasporangiaceae"         "f__Isosphaeraceae"             "f__Legionellaceae"            
#[28] "f__Longimicrobiaceae"          "f__Methylacidiphilaceae"       "f__Methyloligellaceae"        
#[31] "f__Microbacteriaceae"          "f__Micrococcaceae"             "f__Micromonosporaceae"        
#[34] "f__Micropepsaceae"             "f__Microscillaceae"            "f__Moraxellaceae"             
#[37] "f__Myxococcaceae"              "f__Nitrosomonadaceae"          "f__Nitrososphaeraceae"        
#[40] "f__Nocardioidaceae"            "f__Opitutaceae"                "f__Oxalobacteraceae"          
#[43] "f__Parachlamydiaceae"          "f__Pedosphaeraceae"            "f__Phormidiaceae"             
#[46] "f__Phycisphaeraceae"           "f__Pirellulaceae"              "f__Planococcaceae"            
#[49] "f__Pyrinomonadaceae"           "f__Reyranellaceae"             "f__Rhizobiaceae"              
#[52] "f__Rhizobiales_Incertae_Sedis" "f__Rhodanobacteraceae"         "f__Rhodobacteraceae"          
#[55] "f__Rhodospirillaceae"          "f__Rhodothermaceae"            "f__Rickettsiaceae"            
#[58] "f__Roseiflexaceae"             "f__Rubinisphaeraceae"          "f__Saccharimonadaceae"        
#[61] "f__Sandaracinaceae"            "f__Saprospiraceae"             "f__Silvanigrellaceae"         
#[64] "f__Solirubrobacteraceae"       "f__Sphingomonadaceae"          "f__Sporichthyaceae"           
#[67] "f__Steroidobacteraceae"        "f__Tepidisphaeraceae"          "f__Thermoanaerobaculaceae"    
#[70] "f__Vampirovibrionaceae"        "f__Vermiphilaceae"             "f__Verrucomicrobiaceae"       
#[73] "f__Vicinamibacteraceae"        "f__Xanthobacteraceae"          "f__Xanthomonadaceae"          
#[76] "f__Yersiniaceae"        


# Aurantisolimonas, in 	Family:	Chitinophagaceae
# g__Edaphobaculum under family 	Chitinophagaceae
# Heliimonas genus under family Chitinophagaceae
# genus "g__Pseudoflavitalea" to family Chitinophagaceae
# Taibaiella to Chitinophagaceae
k# g__Babeliales")    under phyla = depentiae 

# Genus == "g__Brevifollis" under family   f__Verrucomicrobiaceae
# Genus == "g__Chthonomonadales")  under family  Geobacteraceae
# Genus == "g__Captivus"  under family Paracaedibacteraceae
# Enhydrobacter under Moraxellaceae
# Halocella genus under family 	Halanaerobiaceae
# Genus "g__Ovatusbacter" under Family o__Gammaproteobacteria
# add genus Obscuribacteraceae under c__Vampirivibrionia o__Obscuribacterales f__Obscuribacteraceae
# genus Rhabdanaerobium under 	Eubacteriaceae
# genus Roseisolibacter under Gemmatimonadaceae
# genus Peredibacter to family  Bacteriovoracaceae
# genus Paeniclostridium to family 	Clostridiaceae
# Pedomicrobium to 	Hyphomicrobiaceae
# genus  Phaselicystis which is in  family   Phaselicystidaceae, order Polyangiales
# genus  Tellurimicrobium in family Blastocatellaceae;



##################################### HEATMAP by phylogeny  ###############

otu<-select(otu, -Phyla, -Domain, -Class, -Order, -Family, -Species, -OTU)
n<-c("Genus", "Root4", "Nodule4", "Root1",  "Nodule1",  "Root2",  "Nodule2", "Root3", "Nodule3", "Rhizosphere4", "Rhizosphere1", "Rhizosphere2", "Rhizosphere3") 
colnames(otu)<-n
head(otu)

## filter for taxa that were not present in the plant.
#otu1<-filter(otu, Root4 > 0 | Nodule4 > 0 | Root1 > 0 | Nodule1 > 0 | Root2 > 0 | Nodule2 > 0 | Root3 > 0 | Nodule3 > 0  )

tree.short
head(otu)

#get the otus that are in the the tree.short
otu<-filter(otu, Genus %in% tree.short$tip.label)

# grab correct order 
target<-tree.short$tip.label
otu<-otu[match(target, otu$Genus),]
dim(otu)
# rm Genus column from df
n <- otu$Genus
otu<-subset(otu, select = c( -Genus))
row.names(otu) <- n


############### make matrix
m<-as.matrix(otu)
row.names(m)
m<-structure(m)
# summary stats
max(m)
min(m)
mean(m)
summary(m)
# i think this matrix needs a log transform
m<-log10(m)
m[m== "-Inf"] <- 0
m
summary(m)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/heatmap1.svg", width = 10, height=10 )
# And make the plot with phylogeny
windows(16,12)
heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
dev.off()
# change order
tr<-read.tree(text =  "((Nodule1:3,Nodule2:3, Nodule3:3,Nodule4:3):2, (Root1:3,Root2:3,Root3:3,Root4:3):2, (Rhizosphere1:3, Rhizosphere2:3,Rhizosphere3:3, Rhizosphere4:3):2);")
heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)
dev.off
# no phylogeny
#make a dendrogram              
row_dendro = as.dendrogram(hclust(dist((m))))
windows(10,10)
heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = tr)
dev.off


############################## HEAT MAP group my habitat use type ######

#what do you need?
###genus level data
###df of abundance of active taxa
###habitat use type assigned

dim(otu)
head(otu)
#otu<-select(otu, -Phyla, -Domain, -Class, -Order, -Family, -Species, -OTU)
n<-c("Genus", "Root4", "Nodule4", "Root1",  "Nodule1",  "Root2",  "Nodule2", "Root3", "Nodule3", "Rhizosphere4", "Rhizosphere1", "Rhizosphere2", "Rhizosphere3") 
colnames(otu)<-n
head(otu)

#summarize by location
otu$rhizo <- rowSums(otu %>% dplyr::select(contains("Rhizo"))) 

otu$roots <- rowSums(otu %>% dplyr::select(contains("Root"))) 

otu$nodule <- rowSums(otu %>% dplyr::select(contains("Nod")))

df<-otu%>% select(Genus, rhizo, roots, nodule) %>% filter(rhizo!=0 |nodule!=0 | roots!=0)

head(df)
dim(df)
#df<-t(df)
#df<-as.data.frame(df)
#dim(df)
#df
habitat<-rep(NA, 348)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="rhizo specailist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="root specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="rhizo and nodule generalist"

habitat
unique(habitat)
df$habitat <- habitat
head(df)


# what if we made a little tree with the habitat types
head(df)
t<-df%>% select(Genus, habitat)
head(t)
t<-t[,c(2,1)]

row.names(t)<- NULL
head(t)
## recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}

t
tree<-df2newick(t)
mytree <- read.tree(text=tree)
plot(mytree)
#
############### make matrix
df<-df %>% select(rhizo,roots, nodule)

n<-row.names(df)
#row.names(df)<-habitat

m<-as.matrix(df)
row.names(m)
m<-structure(m)
# summary stats
max(m)
min(m)
mean(m)
summary(m)
# i think this matrix needs a log transform
m<-log10(m)
m[m== "-Inf"] <- 0
m
summary(m)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/heatmap2.svg", width = 10, height=10 )
# And make the plot with phylogeny
windows(16,12)
heatmap.phylo(x = m, Rowp = mytree, Colp = as.phylo(as.hclust(col_dendro)))
dev.off()
# change order
#tr<-read.tree(text =  "((Nodule1:3,Nodule2:3, Nodule3:3,Nodule4:3):2, (Root1:3,Root2:3,Root3:3,Root4:3):2, (Rhizosphere1:3, Rhizosphere2:3,Rhizosphere3:3, Rhizosphere4:3):2);")
heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)
dev.off
# no phylogeny
#make a dendrogram              
row_dendro = as.dendrogram(hclust(dist((m))))
windows(10,10)
heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = tr)
dev.off






############## use the package that determined the phylogenic signal


library(phylosignal)
#install.packages("phylosignal")
library(adephylo)
#library(ape)
library(phylobase)
data(carni19)

#get the tree
tre <- read.tree(text=carni19$tre)
tree.short
#And we create a dataframe of 3 traits for the 19 carnivora species
#. - Body mass - Random values - 
#Simulated values under a Brownian Motion model along the tree
dat <- list()
dat$mass <- carni19$bm
dat$random <- rnorm(19, sd = 10)
dat$bm <- rTraitCont(tre)
dat <- as.data.frame(dat)

#And we create a dataframe of 3 traits for the  153 microbial species
dim(otu)
head(otu)
#rhizosphere - root - nodule - - Random values - -Simulated values under a Brownian Motion model 
##
dat <- list()

dat$rhizo.mean <- rowMeans(otu %>% 
  dplyr::select(contains("Rhizo"))) %>% 
  glimpse()

#dat$rhizo.log <- rowMeans(otu %>% 
#            dplyr::select(contains("Rhizo"))) %>% 
#            log10() %>%
#            glimpse() 

dat$root.mean <- rowMeans(otu %>% 
      dplyr::select(contains("Root"))) %>% 
      glimpse()

#dat$root.log <- rowMeans(otu %>% 
#                dplyr::select(contains("Root"))) %>% 
#                log10() %>%
#                glimpse()

dat$nodule.mean <- rowMeans(otu %>% 
      dplyr::select(contains("Nod"))) %>% 
      glimpse()

#dat$nodule.log <- rowMeans(otu %>% 
#                    dplyr::select(contains("Nod"))) %>% 
#                    log10() %>%
#                    glimpse()


dat$random <- rnorm(153, sd = 10)
dat$bm <- rTraitCont(tree.short)
#turn it into a data frame now?
dat <- as.data.frame(dat)
 
dat[dat== "-Inf"] <- 0 # change -inf into 0
head(dat)

#dim(tree.short)
tree.short
#We can combine phylogeny and traits into a phylo4d object.
p4d <- phylo4d(tree.short, dat)
barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = TRUE)
#
phyloSignal(p4d = p4d, method = "all")
# calculating the phylogenitic signal with all the mehtods

phylosim <- phyloSim(tree = tree.short, method = "all", nsim = 100, reps = 99)
windows(10,10)
plot(phylosim, stacked.methods = FALSE, quantiles = c(0.05, 0.95))
# they are simillar but a little dif
plot.phylosim(phylosim, what = "pval", stacked.methods = TRUE)
# calc the phylogenic signal for the rhizo abundance trait
rhizo.gen <- phyloCorrelogram(p4d, trait = "rhizo.mean")
# root
root.gen <- phyloCorrelogram(p4d, trait = "root.mean")
# nodule
nodule.gen <- phyloCorrelogram(p4d, trait = "nodule.mean")
# for the random trait
window(10,10)
random.gen <- phyloCorrelogram(p4d, trait = "random")
# and fro Simulated values under a Brownian Motion model along the tree
#where continuous traits evolve randomly over time along a branch,
#with a fixed rate. As soon as descents split at a node of the phylogeny,
#evolution on both branches becomes independent
# so there is autocorrelation in a branch
# this a potential null hypthesis other than random
bm.gen <- phyloCorrelogram(p4d, trait = "bm")

plot(rhizo.gen, main= "rhizosphere")
plot(root.gen, main="root")
windows(10,10)
  plot(nodule.gen, main="nodule")

plot(random.gen, main = "random")
plot(bm.gen, main = "Brownian Motion model")


# what is there are certain parts of the taxa this is true for
carni.lipa <- lipaMoran(p4d)
carni.lipa.p4d <- lipaMoran(p4d, as.p4d = TRUE)

barplot.phylo4d(p4d, bar.col=(carni.lipa$p.value < 0.05) + 1, center = FALSE , scale = FALSE)




#####--------heatmap number one just activity famliy level-------########

## select active taxa
ps1<- subset_samples(ps, Fraction=="BONCAT_Active" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
# 1903 taxa from normalizwed by ncells data data
otu<-as.data.frame(t(as.data.frame(otu_table(ps1))))
dim(otu)
# check distrubution

####### agregate to the family level
n<-row.names(otu)
otu<-mutate(otu, OTU=n)
n<-row.names(taxon)
taxon<- mutate(taxon, OTU=n)
otu<-left_join(otu, taxon)
otu<-select(otu, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)
# aggregate some families that have old nomincalture
otu$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu$Family)
otu<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = otu, FUN = sum, na.rm = TRUE)
row.names(otu) <- otu$Family
# rm family col
otu<-select(otu, -Family)
dim(otu)
colnames(otu)
rownames(otu)
n<-c("Root4", "Nodule4", "Root1",  "Nodule1",  "Root2",  "Nodule2", "Root3", "Nodule3", "Rhizosphere4", "Rhizosphere1", "Rhizosphere2", "Rhizosphere3") 
colnames(otu)<-n
head(otu)
## filter for taxa that were not present in the plant.
otu1<-filter(otu, Root4 > 0 | Nodule4 > 0 | Root1 > 0 | Nodule1 > 0 | Root2 > 0 | Nodule2 > 0 | Root3 > 0 | Nodule3 > 0  )
dim(otu1)
head(otu1)
#add back family column
n<-row.names(otu1)
otu1<-mutate(otu1, Family=n)

# grab phylogeny and Read in the tree file 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")
length(tree$tip.label) # look at the tip labels 
# modify tip labels
tree$tip.label <- gsub("\\_1","",tree$tip.label)
tree$tip.label <- gsub("\\'","",tree$tip.label)
#lol my taxa asignment all have a space so I'll remove thAT
otu1$Family <- gsub(" ", "", otu1$Family)
#modeify tip labls
otu1$Family <- gsub('\\_()', '' , otu1$Family)
otu1$Family <- gsub('_\\(*)*', '' , otu1$Family)
otu1$Family <- gsub('*Subgroup1*', '' , otu1$Family)
otu1$Family <- gsub('\\(SR1)', '' , otu1$Family)
otu1$Family <- gsub('*\\(*)*', '' , otu1$Family)
otu1$Family <- gsub ("Abditibacteriaceae"  , "Actinomycetaceae",  otu1$Family ) #in the same clades

length(intersect(unique(otu1$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
length(setdiff(unique(otu1$Family), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
mynames<-(setdiff(unique(otu1$Family), tree$tip.label))
#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, otu1$Family) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE,  cex = .5)
# grab correct order 
target<-tree.short$tip.label
otu1<-otu1[match(target, otu1$Family),]
dim(otu1)
# rm family column from df
n <- otu1$Family
otu1<-subset(otu1, select = c( -Family))
row.names(otu1) <- n


# make matrix
m<-as.matrix(otu1)
row.names(m)
m<-structure(m)
# summary stats
max(m)
min(m)
mean(m)
summary(m)
# i think this matrix needs a log transform
m<-log10(m)
m[m== "-Inf"] <- 0
m
summary(m)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/heatmap1.svg", width = 10, height=10 )
# And make the plot with phylogeny
windows(16,12)
heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
dev.off()
# change order
tr<-read.tree(text =  "((Nodule1:3,Nodule2:3, Nodule3:3,Nodule4:3):2, (Root1:3,Root2:3,Root3:3,Root4:3):2, (Rhizosphere1:3, Rhizosphere2:3,Rhizosphere3:3, Rhizosphere4:3):2);")
heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)
dev.off
# no phylogeny
#make a dendrogram              
row_dendro = as.dendrogram(hclust(dist((m))))
windows(10,10)
heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = tr)
dev.off
        
        



##------- heatmap #2 active plant/ active rhizopshere -------#######
# compare active rhizosphere to active in plant #
# shows what microbes are more active in the plant than the soil

## select active taxa
ps1<- subset_samples(ps, Fraction=="BONCAT_Active" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
# 1903 taxa from normalizwed by ncells data data
# select active taxa in plant
ps.plant<- subset_samples(ps1,Fraction =="BONCAT_Active" & Compartment!="Rhizosphere")
otu.plant<-as.data.frame(t(as.data.frame(otu_table(ps.plant))))
head(otu.plant)
colnames(otu.plant)
# samples C10E, C10N, C1E, C1N, C2E, C2E, C2N, C5E, C5N, C7E, C7N  
# drop sample C7 because it doesnt have a pair in the soil
otu.plant<-subset(otu.plant, select = -c(C7E.POS_S34,C7N.POS_S64))

####### agregate to the family level
n<-row.names(otu.plant)
otu.plant<-mutate(otu.plant, OTU=n)
n<-row.names(taxon)
taxon<- mutate(taxon, OTU=n)
otu.plant<-left_join(otu.plant, taxon)
otu.plant<-select(otu.plant, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)

# aggregate some families that have old nomincalture
otu.plant$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.plant$Family)

otu.plant<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63) ~ Family, data = otu.plant, FUN = sum, na.rm = TRUE)
row.names(otu.plant) <- otu.plant$Family

# rm family col
otu.plant<-select(otu.plant, -Family)
dim(otu.plant)
colnames(otu.plant)
rownames(otu.plant)


## select active taxa in rhizosphere
ps.soil<- subset_samples(ps1,Fraction =="BONCAT_Active" & Compartment=="Rhizosphere")
otu.soil<-as.data.frame(t(as.data.frame(otu_table(ps.soil))))
head(otu.soil)
colnames(otu.soil)

####### agregate to the family level
n<-row.names(otu.soil)
otu.soil<-mutate(otu.soil, OTU=n)
n<-row.names(taxon)
taxon<- mutate(taxon, OTU=n)
otu.soil<-left_join(otu.soil, taxon)
otu.soil<-select(otu.soil, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)
### rm old family name
otu.soil$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.soil$Family)


otu.soil<-aggregate(cbind(C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = otu.soil, FUN = sum, na.rm = TRUE)
dim(otu.soil)
row.names(otu.soil) <- otu.soil$Family

# rm family col
otu.soil<-select(otu.soil, -Family)

# samples C10R, C1R, C2R, C5R
# make df of active in rhizosphere to divid by
df<-cbind(otu.soil, otu.soil)
df<-df[sort(colnames(df))]
colnames(df)
dim(df)
dim(otu.plant)
#hist(df1$C10E.POS_S60, breaks = 20)
# add 1 to everything (absent in active & absent in total = no change)
otu_log2<-log2(otu.plant/(df+.1))

dim(otu_log2)
colnames(otu_log2)
n<-c("C10E", "C10N","C1E" , "C1N" , "C2E" , "C2N" , "C5E" , "C5N")
colnames(otu_log2)<-n
head(otu_log2)
# check distribution
otu_log2
#hist(otu_log2$C10N) # most things didn't change because most taxa not present 
#hist(otu_log2$C5N)
hist(otu_log2$C1E, breaks = 20)
# col dendogram doesn't work because a ton of stuff equals 0 infinit

## filter for taxa that were not present in the plant.

otu_log2[otu_log2== "-Inf"] <- -11

otu_log2.1<-filter(otu_log2, C10E > -11 | C10N > -11 | C1E > -11 | C1N > -11 | C2E > -11 | C2N > -11 | C5E > -11 | C5N > -11  )

dim(otu_log2.1)
otu_log2.1

# we re code -inf as -99
# this means a taxa wasn't present in the location sampled

# build your data your trying to pu tin the tree
n<-row.names(otu_log2.1)
otu_log2.1<-mutate(otu_log2.1, Family=n)

############## make heatmap #2-------#################

# grab phylogeny and Read in the tree file 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")


length(tree$tip.label) # look at the tip labels 

#this is how to modify tip labels
tree$tip.label <- gsub("\\_1","",tree$tip.label)
tree$tip.label <- gsub("\\'","",tree$tip.label)
#lol my taxa asignment all have a space so I'll remove thAT
otu_log2.1$Family <- gsub(" ", "", otu_log2.1$Family)
#modeify tip labls
otu_log2.1$Family <- gsub('\\_()', '' , otu_log2.1$Family)
otu_log2.1$Family <- gsub('_\\(*)*', '' , otu_log2.1$Family)
otu_log2.1$Family <- gsub('*Subgroup1*', '' , otu_log2.1$Family)
otu_log2.1$Family <- gsub('\\(SR1)', '' , otu_log2.1$Family)
otu_log2.1$Family <- gsub('*\\(*)*', '' , otu_log2.1$Family)
otu_log2.1$Family <- gsub ("Abditibacteriaceae"  , "Actinomycetaceae",  otu_log2.1$Family ) #in the same clades


length(intersect(unique(otu_log2.1$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
length(setdiff(unique(otu_log2.1$Family), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
mynames<-(setdiff(unique(otu_log2.1$Family), tree$tip.label))

# CHECK NAMES IN BIG TREE
#spnames = tree$tip.label
#spnames<-sort(spnames)
#spnames[grep("^[P][r].*", spnames)]

# CHECK NAEMS OF TAXA THAT ARE MISSIG
#mynames<-sort(mynames)
#mynames
#mynames[grep("^[B].*", mynames)]
#out_log2.1$Family[grep("^[A].*", out_log2.1$Family)]
##### MAKE TREE SMALLER

# shorten phylogeny to match what is in our log 2 file
#families we want
asvs_remove<-setdiff(tree$tip.label, otu_log2.1$Family) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need

plot(tree.short, no.margin=TRUE,  cex = .5)
#nodelabels()
#### add tip "Blastocatellaceae"  to "Acidobacteriaceae" 
#tree.short<-AddTip(tree.short, where = 118, "Blastocatellaceae")
### phylum Latescibacterota 
### Prevotellaceae family in order Bacteriodales


# order the table so it matches the 
# make otu column 
# grab correct order 
target<-tree.short$tip.label
head(otu_log2.1)
dim(otu_log2.1)
otu_log2.1<-otu_log2.1[match(target, otu_log2.1$Family),]
dim(otu_log2.1)
# rm out column
# makes sure no space in col names
n <- otu_log2.1$Family
otu_log2.1<-subset(otu_log2.1, select = c( -Family))
row.names(otu_log2.1) <- n
                 
colnames(otu_log2.1)
head(otu_log2.1)
dim(otu_log2.1)
tree.short
#row.names(out_log2.1) <-f

# make matrix
     m<-as.matrix(otu_log2.1)
      row.names(m)
      #row.names(m)<-f # make row names family names
      # Create the matrix and get the column dendrogram for the heatmap from it.
      m<-structure(m)
      #make a dendrogram              
      col_dendro = as.dendrogram(hclust(dist(t(m))))
      # 
      n <- colnames(m) 
      
 # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(14,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp =as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
         Rowp = tree.short
         Colp = as.phylo(as.hclust(col_dendro))
         #heatmap.phylo <- function(x, Rowp, Colp, ...) {
          # l = length(seq(-4.9, 5, 0.1))
          # pal = colorRampPalette(c('#2166ac', '#f7f7f7', '#b2182b'))(l)
          
  #make a plot without phylogeny 
        #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(12,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         

#____________________________________________#
################# active location X /totallocation X  heatmap #3############
#### value in active compared to total for that respective location

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
   
  ####### agregate to the family level
  n<-row.names(otu.total)
  otu.total<-mutate(otu.total, OTU=n)
  n<-row.names(taxon)
  taxon<- mutate(taxon, OTU=n)
  otu.total<-left_join(otu.total, taxon)
  otu.total<-select(otu.total, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)
  
  # aggregate some families that have old nomincalture
  otu.total$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.total$Family)
  
  otu.total<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13  ,C1R.SYBR_S16,  C2E.SYBR_S22 ,
                             C2N.SYBR_S15,  C2R.SYBR_S17,  C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24 , C7N.SYBR_S25 ) ~ Family, data = otu.total, FUN = sum, na.rm = TRUE)
  row.names(otu.total) <- otu.total$Family
  
  # rm family col
  otu.total<-select(otu.total, -Family)
  dim(otu.total)
  colnames(otu.total)
  rownames(otu.total)
  

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
    
    ####### agregate to the family level
    n<-row.names(otu.active)
    otu.active<-mutate(otu.active, OTU=n)
    n<-row.names(taxon)
    taxon<- mutate(taxon, OTU=n)
    otu.active<-left_join(otu.active, taxon)
    otu.active<-select(otu.active, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)
    
    # aggregate some families that have old nomincalture
    otu.active$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.active$Family)
    
    otu.active<-aggregate(cbind(C10N.POS_S65, C10R.POS_S30, C1E.POS_S31,  C1N.POS_S61,  C1R.POS_S27,  C2E.POS_S32, 
                                 C2N.POS_S62,  C2R.POS_S28, C5E.POS_S33,  C5R.POS_S29,  C7E.POS_S34,  C7N.POS_S64  ) ~ Family, data = otu.active, FUN = sum, na.rm = TRUE)
    row.names(otu.active) <- otu.active$Family
    
    # rm family col
    otu.active<-select(otu.active, -Family)
    dim(otu.active)
    colnames(otu.active)
    rownames(otu.active)    
   
 # log fold change from active to inactive 
 # add 1 to everything absent in total = no change)
         
         otu_log2<-log2(otu.active/(otu.total+.1))
         
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
         # we re code -inf as -11
         # this means a taxa wasn't present in the location sampled
         otu_log2[otu_log2== "-Inf"] <- -11
           ### remove taxa not present in plant
          colnames(otu_log2)
         #otu_log2<- filter(otu_log2, C10N > -11 | C1E > -11 | C1N > -11 | C2E > -11 | C2N > -11 | C5E > -11 | C7E > -11 | C7N > -11  )
          # add back family column
         f<-row.names(otu_log2)
         otu_log2<-mutate(otu_log2, Family = f)
         # remove columns that are not rhizo for heat map #4
         #otu_log2.r<-otu_log2 %>% select(-C10N  , -C1E  , -C1N  , -C2E  , -C2N  , -C5E  , -C7E  , -C7N)
         
########### heatmap #3 make the tree -------#################
         # grab phylogeny and Read in the tree file 
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
         tree = read.tree("family.nwk")
         
         # chec if the ncbi one is better
         #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/NCBI")
         #tree = read.tree("genus.nwk")
         # it is not 
         length(tree$tip.label) # look at the tip labels 
         
         #this is how to modify tip labels
         tree$tip.label <- gsub("\\_1","",tree$tip.label)
         tree$tip.label <- gsub("\\'","",tree$tip.label)
         #lol my taxa asignment all have a space so I'll remove thAT
         otu_log2$Family <- gsub(" ", "", otu_log2$Family)
         #modeify tip labls
         otu_log2$Family <- gsub('\\_()', '' , otu_log2$Family)
         otu_log2$Family <- gsub('_\\(*)*', '' , otu_log2$Family)
         otu_log2$Family <- gsub('*Subgroup1*', '' , otu_log2$Family)
         otu_log2$Family <- gsub('\\(SR1)', '' , otu_log2$Family)
         otu_log2$Family <- gsub('*\\(*)*', '' , otu_log2$Family)
         otu_log2$Family <- gsub ("Abditibacteriaceae"  , "Actinomycetaceae",  otu_log2$Family ) #in the same clades
         
         
         length(intersect(unique(otu_log2$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
         length(setdiff(unique(otu_log2$Family), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
         mynames<-(setdiff(unique(otu_log2$Family), tree$tip.label))
         
         # CHECK NAMES IN BIG TREE
         #spnames = tree$tip.label
         #spnames<-sort(spnames)
         #spnames[grep("^[P][r].*", spnames)]
         
         # CHECK NAEMS OF TAXA THAT ARE MISSIG
         mynames<-sort(mynames)
         mynames
         mynames[grep("^[Aa].*", mynames)]
         otu_log2$Family[grep("^[T].*", otu_log2$Family)]
         ##### MAKE TREE SMALLER
         
         # shorten phylogeny to match what is in our log 2 file
         #families we want
         asvs_remove<-setdiff(tree$tip.label, otu_log2$Family) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
        # windows(14,14)
         #plot(tree.short, no.margin=TRUE,  cex = .5)
         #nodelabels()
         #ape::nodelabels( 253, 253, bg="green")
         #ape::nodelabels(31, 31, bg="green")
         
         ### rename node labels
         tree.short[["node.label"]] <- c(1:138)
         
         # add branches    
         #which(tree.short$tip.label=="Acidobacteriaceae")
         #which(tree.short$tip.label== "Blastocatellia"  )
         
         #88
        # tree.short$node.label[88]
         #tree.short<-AddTip(tree.short, where = 88, "Blastocatellia"  )
         #ree.short<-AddTip(tree.short, where = 146,  "Blastocatellaceae" )
         # 146
         #### add branches
         #which(tree.short$tip.label=="Acidobacteriaceae")
         ##
         # add up higher
         #tree.short<-AddTip(tree.short, where = 234,  "Acidobacteriales") # order
         #tree.short<-AddTip(tree.short, where = 234, "Acidobacteriae" )  #class?
         
         # add branches    
         #which(tree.short$tip.label== "Rhizobiaceae")
         # 122
         ## branch at node 251
         #tree.short<-AddTip(tree.short, where = 253,  "Alphaproteobacteria")
         #which(tree.short$tip.label== "Alphaproteobacteria" )
         #### add branch "Actinobacteria" )
         #which(tree.short$tip.label=="Actinomycetaceae" )
         # MAYBE 158
         #tree.short<-AddTip(tree.short, where = 158, "Actinobacteria"   )
         
         
         #match with short tree
         length(intersect(unique(otu_log2$Family), tree.short$tip.label)) # Apply setdiff function to see what's missing from the tree
         length(setdiff(unique(otu_log2$Family), tree.short$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
         mynames<-(setdiff(unique(otu_log2$Family), tree.short$tip.label))
         
         
         
         ### phylum Latescibacterota 
         ### Prevotellaceae family in order Bacteriodales
         
         
         
         # order the table so it matches the 
         # make otu column 
         # grab correct order 
         target<-tree.short$tip.label
         head(otu_log2)
         dim(otu_log2)
         otu_log2<-otu_log2[match(target, otu_log2$Family),]
         dim(otu_log2)
         # rm out column
         # makes sure no space in col names
         n <- otu_log2$Family
         otu_log2<-subset(otu_log2, select = c( -Family))
         row.names(otu_log2) <- n
         
         colnames(otu_log2)
         head(otu_log2)
         dim(otu_log2)
         tree.short
         #row.names(out_log2.1) <-f
         
        
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
         windows(12,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))

         
         
         ######### change order of samples in tree
        as.phylo(colnames(m))
        Colp<-read.tree(text =  "((C10N:3,C1N:3):2, (C1E:3,C1N:3,C1R:3):2,   (C2E:3,C2N:3,C2R:3):2, (C5E:3,C5R:3):2, (C7E:3  ,C7N:3):2);")

        

        # And make the plot with phylogeny
        #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
        windows(12,10)
        heatmap.phylo(x = m, Rowp = tree.short, Colp = Colp)
        #dev.off()
         
###### subset taxa present in plant and change order #############
        
        otu_log2
        otu_log2<-subset(otu_log2, otu_log2$C1E != -11 | otu_log2$C2E != -11  | otu_log2$C5E != -11 | otu_log2$C7E != -11 )
        f<-row.names(otu_log2)
        otu_log2<-mutate(otu_log2, Family = f)
        
        
########### heatmap #3 make the tree -------#################
        # grab phylogeny and Read in the tree file 
        setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
        tree = read.tree("family.nwk")
        
        # chec if the ncbi one is better
        #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/NCBI")
        #tree = read.tree("genus.nwk")
        # it is not 
        length(tree$tip.label) # look at the tip labels 
        
        #this is how to modify tip labels
        tree$tip.label <- gsub("\\_1","",tree$tip.label)
        tree$tip.label <- gsub("\\'","",tree$tip.label)
        #lol my taxa asignment all have a space so I'll remove thAT
        otu_log2$Family <- gsub(" ", "", otu_log2$Family)
        #modeify tip labls
        otu_log2$Family <- gsub('\\_()', '' , otu_log2$Family)
        otu_log2$Family <- gsub('_\\(*)*', '' , otu_log2$Family)
        otu_log2$Family <- gsub('*Subgroup1*', '' , otu_log2$Family)
        otu_log2$Family <- gsub('\\(SR1)', '' , otu_log2$Family)
        otu_log2$Family <- gsub('*\\(*)*', '' , otu_log2$Family)
        otu_log2$Family <- gsub ("Abditibacteriaceae"  , "Actinomycetaceae",  otu_log2$Family ) #in the same clades
        
        
        length(intersect(unique(otu_log2$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
        length(setdiff(unique(otu_log2$Family), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
        mynames<-(setdiff(unique(otu_log2$Family), tree$tip.label))
        
    
        
        # CHECK NAEMS OF TAXA THAT ARE MISSIG
        mynames<-sort(mynames)
        mynames
        mynames[grep("^[Aa].*", mynames)]
        otu_log2$Family[grep("^[T].*", otu_log2$Family)]
        ##### MAKE TREE SMALLER
        
        # shorten phylogeny to match what is in our log 2 file
        #families we want
        asvs_remove<-setdiff(tree$tip.label, otu_log2$Family) #asvs we don't want
        tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
        # windows(14,14)
      
        ### rename node labels
        tree.short[["node.label"]] <- c(1:138)
        
        #match with short tree
        length(intersect(unique(otu_log2$Family), tree.short$tip.label)) # Apply setdiff function to see what's missing from the tree
        length(setdiff(unique(otu_log2$Family), tree.short$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
        mynames<-(setdiff(unique(otu_log2$Family), tree.short$tip.label))
        
               # order the table so it matches the 
        # make otu column 
        # grab correct order 
        target<-tree.short$tip.label
        head(otu_log2)
        dim(otu_log2)
        otu_log2<-otu_log2[match(target, otu_log2$Family),]
        dim(otu_log2)
        # rm out column
        # makes sure no space in col names
        n <- otu_log2$Family
        otu_log2<-subset(otu_log2, select = c( -Family))
        row.names(otu_log2) <- n
        
        colnames(otu_log2)
        head(otu_log2)
        dim(otu_log2)
        tree.short
        #row.names(out_log2.1) <-f
        
        
        # make matrix
        m<-as.matrix(otu_log2)
        row.names(m)
        #row.names(m)<-f # make row names family names
        # Create the matrix and get the column dendrogram for the heatmap from it.
        m<-structure(m)
        #make a dendrogram              
        col_dendro = as.dendrogram(hclust(dist(t(m))))
        
        
        setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
        
        # And make the plot with phylogeny
        svg(file="figures/16s/heatmap3.svg",width = 10, height=7 )
        
        windows(16,10)
        #heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
        #dev.off()
        
        tr<-read.tree(text =  "((C1N:3,C2N:3,C7N:3,C10N:3):2, (C1E:3,C2E:3,C5E:3,C7E:3):2,   (C1R:3, C2R:3,C5R:3, C10R:3):2);")
       
        heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)
        
        dev.off
        
          heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)
        
        # no phylogeny
        
        #make a plot without phylogeny 
        #make a dendrogram              
        row_dendro = as.dendrogram(hclust(dist((m))))
        tr<-read.tree(text =  "((C1N:3,C2N:3,C7N:3,C10N:3):2, (C1E:3,C2E:3,C5E:3,C7E:3):2,   (C1R:3, C2R:3,C5R:3, C10R:3):2);")

        ############  make the heatmap  plot #################33
        windows(16,10)
        heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = tr)
        
        
        dev.off
        
        
        
        
###subset by just rhizosphere
         colnames(otu_log2)
         df<-subset(otu_log2, select = c(C10R, C1R, C2R, C5R ))
         df<-subset(df, df$C10R != -99 | df$C1R != -99  | df$C2R != -99 | df$C5R != -99 )
         dim(df)
         n<-row.names(df)
         df<-df%>% mutate( Family = n)
         # grab tree
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
         tree = read.tree("family.nwk")
         length(tree$tip.label) # look at the tip labels
         #modify tip labels
         tree$tip.label <- gsub("\\_1","",tree$tip.label)
         tree$tip.label <- gsub("\\'","",tree$tip.label)
         # shorten tree
         asvs_remove<-setdiff(tree$tip.label, df$Family) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         # grab correct order 
         target<-tree.short$tip.label
         head(df)
         df<-df[match(target, df$Family),]
         dim(df)
         # rm out column
         # makes sure no space in col names
         n <- df$Family
         df<-subset(df, select = c( -Family))
         row.names(df) <- n

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
         
#######heat map #3 just taxa in the plant #######
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
         
         Colp<-read.tree(text =  "((C10N:3,C1N:3):2, (C1E:3,C1N:3,C1R:3):2,   (C2E:3,C2N:3,C2R:3):2, (C5E:3,C5R:3):2, (C7E:3  ,C7N:3):2);")
         
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
         
########## heatmaps # 3 active/ total at order level ####################
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
         ####### agregate to the family level
         n<-row.names(otu.total)
         otu.total<-mutate(otu.total, OTU=n)
         n<-row.names(taxon)
         taxon<- mutate(taxon, OTU=n)
         otu.total<-left_join(otu.total, taxon)
         otu.total<-select(otu.total, -Phyla, -Domain, -Class, -Family, -Genus, -Species, -OTU)
         
         # aggregate some families that have old nomincalture
         otu.total$Order <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.total$Order)
         
         otu.total<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13  ,C1R.SYBR_S16,  C2E.SYBR_S22 ,
                                    C2N.SYBR_S15,  C2R.SYBR_S17,  C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24 , C7N.SYBR_S25 ) ~ Order, data = otu.total, FUN = sum, na.rm = TRUE)
         row.names(otu.total) <- otu.total$Order
         
         # rm order col
         otu.total<-select(otu.total, -Order)
         dim(otu.total)
         colnames(otu.total)
         rownames(otu.total)
         
         
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
         
         ####### agregate to the order level
         n<-row.names(otu.active)
         otu.active<-mutate(otu.active, OTU=n)
         n<-row.names(taxon)
         taxon<- mutate(taxon, OTU=n)
         otu.active<-left_join(otu.active, taxon)
         otu.active<-select(otu.active, -Phyla, -Domain, -Class, -Family, -Genus, -Species, -OTU)
         
         # aggregate some families that have old nomincalture
         otu.active$Order <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.active$Order)
         
         otu.active<-aggregate(cbind(C10N.POS_S65, C10R.POS_S30, C1E.POS_S31,  C1N.POS_S61,  C1R.POS_S27,  C2E.POS_S32, 
                                     C2N.POS_S62,  C2R.POS_S28, C5E.POS_S33,  C5R.POS_S29,  C7E.POS_S34,  C7N.POS_S64  ) ~ Order, data = otu.active, FUN = sum, na.rm = TRUE)
         row.names(otu.active) <- otu.active$Order
         
         # rm Order col
         otu.active<-select(otu.active, -Order)
         dim(otu.active)
         colnames(otu.active)
         rownames(otu.active)    
         
         # log fold change from active to inactive 
         # add 1 to everything absent in total = no change)
         
         otu_log2<-log2(otu.active/(otu.total+.1))
         
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
         # we re code -inf as -11
         # this means a taxa wasn't present in the location sampled
         otu_log2[otu_log2== "-Inf"] <- -11
         ### remove taxa not present in plant
         colnames(otu_log2)
         #otu_log2<- filter(otu_log2, C10N > -11 | C1E > -11 | C1N > -11 | C2E > -11 | C2N > -11 | C5E > -11 | C7E > -11 | C7N > -11  )
         # add back family column
         f<-row.names(otu_log2)
         otu_log2<-mutate(otu_log2, Order = f)
         # remove columns that are not rhizo for heat map #4

         ########### heatmap #3 make the tree -------#################
         # grab phylogeny and Read in the tree file 
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
         tree = read.tree("Order.nwk")
         
         # chec if the ncbi one is better
         #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/NCBI")
         #tree = read.tree("genus.nwk")
         # it is not 
         length(tree$tip.label) # look at the tip labels 
         
         #this is how to modify tip labels
         #lol my taxa asignment all have a space so I'll remove thAT
         otu_log2$Order <- gsub(" ", "", otu_log2$Order)
         #modeify tip labls
         length(intersect(unique(otu_log2$Order), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
         length(setdiff(unique(otu_log2$Order), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
         mynames<-(setdiff(unique(otu_log2$Order), tree$tip.label))
         
         # CHECK NAMES IN BIG TREE
         #spnames = tree$tip.label
         #spnames<-sort(spnames)
         #spnames[grep("^[P][r].*", spnames)]
         
         # CHECK NAEMS OF TAXA THAT ARE MISSIG
         mynames<-sort(mynames)
         mynames
         mynames[grep("^[Aa].*", mynames)]
         tree$tip.label[grep("^[Aa].*", tree$tip.label)]
         otu_log2$Order[grep("^[T].*", otu_log2$Order)]
         ##### MAKE TREE SMALLER
         
         # shorten phylogeny to match what is in our log 2 file
         #families we want
         asvs_remove<-setdiff(tree$tip.label, otu_log2$Order) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         windows(14,14)
         plot(tree.short, no.margin=TRUE,  cex = .5)
         nodelabels()
         ape::nodelabels( 253, 253, bg="green")
         ape::nodelabels(31, 31, bg="green")
         
         ### rename node labels
         #tree.short[["node.label"]] <- c(1:138)
         
         # add branches    
         #which(tree.short$tip.label=="Acidobacteriaceae")
         #which(tree.short$tip.label== "Blastocatellia"  )
         
         #88
         # tree.short$node.label[88]
         #tree.short<-AddTip(tree.short, where = 88, "Blastocatellia"  )
         #ree.short<-AddTip(tree.short, where = 146,  "Blastocatellaceae" )
         # 146
         #### add branches
         #which(tree.short$tip.label=="Acidobacteriaceae")
         ##
         # add up higher
         #tree.short<-AddTip(tree.short, where = 234,  "Acidobacteriales") # order
         #tree.short<-AddTip(tree.short, where = 234, "Acidobacteriae" )  #class?
         
         # add branches    
         #which(tree.short$tip.label== "Rhizobiaceae")
         # 122
         ## branch at node 251
         #tree.short<-AddTip(tree.short, where = 253,  "Alphaproteobacteria")
         #which(tree.short$tip.label== "Alphaproteobacteria" )
         #### add branch "Actinobacteria" )
         #which(tree.short$tip.label=="Actinomycetaceae" )
         # MAYBE 158
         #tree.short<-AddTip(tree.short, where = 158, "Actinobacteria"   )
         
         
         #match with short tree
         length(intersect(unique(otu_log2$Order), tree.short$tip.label)) # Apply setdiff function to see what's missing from the tree
         length(setdiff(unique(otu_log2$Order), tree.short$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
         mynames<-(setdiff(unique(otu_log2$Order), tree.short$tip.label))
         
         
         
         ### phylum Latescibacterota 
         ### Prevotellaceae family in order Bacteriodales
         
         
         
         # order the table so it matches the 
         # make otu column 
         # grab correct order 
         target<-tree.short$tip.label
         head(otu_log2)
         dim(otu_log2)
         otu_log2<-otu_log2[match(target, otu_log2$Order),]
         dim(otu_log2)
         # rm out column
         # makes sure no space in col names
         n <- otu_log2$Order
         otu_log2<-subset(otu_log2, select = c( -Order))
         row.names(otu_log2) <- n
         
         colnames(otu_log2)
         head(otu_log2)
         dim(otu_log2)
         tree.short
         #row.names(out_log2.1) <-f
         
         
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
         windows(12,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(8,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
         
         
         ######### change order of samples in tree
         as.phylo(colnames(m))
         Colp<-read.tree(text =  "((C10N:3,C10R:3):2, (C1E:3,C1N:3,C1R:3):2,   (C2E:3,C2N:3,C2R:3):2, (C5E:3,C5R:3):2, (C7E:3  ,C7N:3):2);")
         
         # And make the plot with phylogeny
         #pdf(file="../figures/heatmap.pdf",width = 9,height=10, useDingbats=FALSE)
         windows(12,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = Colp)
         #dev.off()
         
         
         
         
         ###subset by just rhizosphere
         colnames(otu_log2)
         df<-subset(otu_log2, select = c(C10R, C1R, C2R, C5R ))
         df<-subset(df, df$C10R != -99 | df$C1R != -99  | df$C2R != -99 | df$C5R != -99 )
         dim(df)
         n<-row.names(df)
         df<-df%>% mutate( Family = n)
         # grab tree
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
         tree = read.tree("family.nwk")
         length(tree$tip.label) # look at the tip labels
         #modify tip labels
         tree$tip.label <- gsub("\\_1","",tree$tip.label)
         tree$tip.label <- gsub("\\'","",tree$tip.label)
         # shorten tree
         asvs_remove<-setdiff(tree$tip.label, df$Family) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         # grab correct order 
         target<-tree.short$tip.label
         head(df)
         df<-df[match(target, df$Family),]
         dim(df)
         # rm out column
         # makes sure no space in col names
         n <- df$Family
         df<-subset(df, select = c( -Family))
         row.names(df) <- n
         
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
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
                  
#####heat map # 1  total/total ######
# what is the selectivity of the plant?
# total viable in fraction X / total viable in rhizosphere
         
## select total taxa
         ps1<- subset_samples(ps, Fraction=="Total_Cells" ) 
         sample_data(ps)
         ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
         any(taxa_sums(ps1) == 0)
         
         
         ps1
         # 2473 taxa from non rareified data
         #filter for taxa that have at least 10 cells
         ps1<-prune_taxa(taxa_sums(ps1) > 10, ps1)
         ps1
         # 2269 taxa
         
# select total viable cells in plant
         ps.plant<- subset_samples(ps1,Fraction =="Total_Cells" & Compartment!="Rhizosphere")
         otu.plant<-as.data.frame(t(as.data.frame(otu_table(ps.plant))))
         head(otu.plant)
         colnames(otu.plant)
         dim(otu.plant)
         
         
         
         # 8 columns
    
        ####### agregate to the family level
         n<-row.names(otu.plant)
         otu.plant<-mutate(otu.plant, OTU=n)
         n<-row.names(taxon)
         taxon<- mutate(taxon, OTU=n)
         otu.plant<-left_join(otu.plant, taxon)
         
         # if out plant family col in blank then sub for order
         #For each row in data file
         #if otu.plant$Family == "",
         #then value = otu$order
         otu.plant[is.na(otu.plant)] <- ""

         otu.plant <- otu.plant %>% 
           mutate(Family = ifelse(Family == "", Order,Family))
         otu.plant <- otu.plant %>% 
           mutate(Family = ifelse(Family == "", Class,Family))
         
         
        otu.plant<-select(otu.plant, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)
         
         # aggregate some families that have old nomincalture
         otu.plant$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.plant$Family)
         
         colnames(otu.plant)
         otu.plant<-aggregate(cbind(C10N.SYBR_S26, C1E.SYBR_S21,  C1N.SYBR_S13,  C2E.SYBR_S22, C2N.SYBR_S15,  C5E.SYBR_S23,  C7E.SYBR_S24 
                                     , C7N.SYBR_S25 ) ~ Family, data = otu.plant, FUN = sum, na.rm = TRUE)
         row.names(otu.plant) <- otu.plant$Family
         
         # rm family col
         otu.plant<-select(otu.plant, -Family)
         dim(otu.plant)
         # 308 x 8
         colnames(otu.plant)
         rownames(otu.plant)
                  
  
 
## select total taxa in rhizosphere
         ps.soil<- subset_samples(ps1,  Compartment=="Rhizosphere" & Fraction=="Total_Cells" )
         otu.soil<-as.data.frame(t(as.data.frame(otu_table(ps.soil))))
         head(otu.soil)
         colnames(otu.soil)
         dim(otu.soil)
         # samples C10R, C1R, C2R, C5R, C7R
         
         ####### agregate to the family level
         n<-row.names(otu.soil)
         otu.soil<-mutate(otu.soil, OTU=n)
         n<-row.names(taxon)
         taxon<- mutate(taxon, OTU=n)
         otu.soil<-left_join(otu.soil, taxon)
         
         otu.soil[is.na(otu.soil)] <- ""
         
         otu.soil <- otu.soil %>% 
           mutate(Family = ifelse(Family == "", Order,Family))
         otu.soil <- otu.soil %>% 
           mutate(Family = ifelse(Family == "", Class,Family))
         
         otu.soil<-select(otu.soil, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -OTU)
         
         # aggregate some families that have old nomincalture
         otu.soil$Family <- gsub("*Rhizobiales_Incertae_Sedis*", "Rhizobiaceae" , otu.soil$Family)
         
         colnames(otu.soil)
         otu.soil<-aggregate(cbind(C10R.SYBR_S20, C1R.SYBR_S16, C2R.SYBR_S17, C5R.SYBR_S18, C7R.SYBR_S19) ~ Family, data = otu.soil, FUN = sum, na.rm = TRUE)
         row.names(otu.soil) <- otu.soil$Family
         
         # rm family col
         otu.soil<-select(otu.soil, -Family)
         dim(otu.soil)
         # 299 x 5
         colnames(otu.soil)
         rownames(otu.soil)
         
        
         # make otu soil df the roght size to make with the plant df
         df<-cbind(otu.soil, otu.soil)
         colnames(df)
         dim(df)
         # remove 2 columns
         df<-df[, order(colnames(df))]
         df<-subset(df, select= c(-C10R.SYBR_S20.1, -C5R.SYBR_S18.1))
        dim(df)
         
# calculate log fold change         
# add 1 to everything (absent in active & absent in total = no change)
         
         otu_log2<-log2(otu.plant/(df+.1))
         
         dim(otu_log2)
         colnames(otu_log2)
         n<-c("C10N", "C1E","C1N", "C2E" , "C2N" , "C5E" , "C7E", "C7N")
         colnames(otu_log2)<-n
         head(otu_log2)
         # check distribution
         otu_log2
         hist(otu_log2$C10N) # most things didn't change because most taxa not present 
         hist(otu_log2$C7N)
         hist(otu_log2$C1E, breaks = 20)
         # col dendogram doesn't work because a ton of stuff equals 0 infinit
         # we re code -inf as -99
         # this means a taxa wasn't present in the location sampled
         
         
         otu_log2[otu_log2== "-Inf"] <- -15
         ## add family column
         f<-row.names(otu_log2)
         otu_log2<-mutate(otu_log2, Family = f)
         
         #what are the values for these missing taxa
         
         otu_log2<- filter(otu_log2, C10N > -15 | C1E > -15 | C1N > -15 | C2E > -15 | C2N > -15 | C5E > -15 | C7E > -15 | C7N > -15  )
         

         
####### make heatmap # 1 ############
         # grab phylogeny and Read in the tree file 
         
         # grab phylogeny and Read in the tree file 
         setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
         tree = read.tree("family.nwk")
         
         # chec if the ncbi one is better
         #setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/NCBI")
         #tree = read.tree("genus.nwk")
         # it is not 
         length(tree$tip.label) # look at the tip labels 
         
         #this is how to modify tip labels
         tree$tip.label <- gsub("\\_1","",tree$tip.label)
         tree$tip.label <- gsub("\\'","",tree$tip.label)
         #lol my taxa asignment all have a space so I'll remove thAT
         otu_log2$Family <- gsub(" ", "", otu_log2$Family)
         #modeify tip labls
         otu_log2$Family <- gsub('\\_()', '' , otu_log2$Family)
         otu_log2$Family <- gsub('_\\(*)*', '' , otu_log2$Family)
         otu_log2$Family <- gsub('*Subgroup1*', '' , otu_log2$Family)
         otu_log2$Family <- gsub('\\(SR1)', '' , otu_log2$Family)
         otu_log2$Family <- gsub('*\\(*)*', '' , otu_log2$Family)
         otu_log2$Family <- gsub ("Abditibacteriaceae"  , "Actinomycetaceae",  otu_log2$Family ) #in the same clades
         
         
         length(intersect(unique(otu_log2$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
         length(setdiff(unique(otu_log2$Family), tree$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
         mynames<-(setdiff(unique(otu_log2$Family), tree$tip.label))
         
         # CHECK NAMES IN BIG TREE
         #spnames = tree$tip.label
         #spnames<-sort(spnames)
         #spnames[grep("^[P][r].*", spnames)]
         
         # CHECK NAEMS OF TAXA THAT ARE MISSIG
         mynames<-sort(mynames)
         mynames
         mynames[grep("^[Aa].*", mynames)]
         otu_log2$Family[grep("^[T].*", otu_log2$Family)]
         ##### MAKE TREE SMALLER
         #what are the values for these missing taxa 
         df1<-otu_log2[otu_log2$Family %in% mynames,] %>%
          filter( C10N > -11 | C1E > -11 | C1N > -11 | C2E > -11 | C2N > -11 | C5E > -11 | C7E > -11 | C7N > -11  )
         
         mynames<-df1$Family
         
       # shorten phylogeny to match what is in our log 2 file
         #families we want
         asvs_remove<-setdiff(tree$tip.label, otu_log2$Family) #asvs we don't want
         tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
         windows(14,14)
         plot(tree.short, no.margin=TRUE,  cex = .5)
          nodelabels()
          ape::nodelabels( 253, 253, bg="green")
          ape::nodelabels(31, 31, bg="green")
          
        ### rename node labels
        tree.short[["node.label"]] <- c(1:138)

      # add branches    
        which(tree.short$tip.label=="Acidobacteriaceae")
        which(tree.short$tip.label== "Blastocatellia"  )
        
        #88
        tree.short$node.label[88]
        tree.short<-AddTip(tree.short, where = 88, "Blastocatellia"  )
        tree.short<-AddTip(tree.short, where = 146,  "Blastocatellaceae" )
        # 146
        #### add branches
        which(tree.short$tip.label=="Acidobacteriaceae")
        ##
        # add up higher
        tree.short<-AddTip(tree.short, where = 234,  "Acidobacteriales") # order
        tree.short<-AddTip(tree.short, where = 234, "Acidobacteriae" )  #class?

        # add branches    
        which(tree.short$tip.label== "Rhizobiaceae")
        # 122
        ## branch at node 251
        tree.short<-AddTip(tree.short, where = 253,  "Alphaproteobacteria")
        which(tree.short$tip.label== "Alphaproteobacteria" )
        #### add branch "Actinobacteria" )
        which(tree.short$tip.label=="Actinomycetaceae" )
        # MAYBE 158
        tree.short<-AddTip(tree.short, where = 158, "Actinobacteria"   )
        

        #match with short tree
        length(intersect(unique(otu_log2$Family), tree.short$tip.label)) # Apply setdiff function to see what's missing from the tree
        length(setdiff(unique(otu_log2$Family), tree.short$tip.label))# Apply setdiff function to see what's missing from tree in but in my df
        mynames<-(setdiff(unique(otu_log2$Family), tree.short$tip.label))
        
        

         ### phylum Latescibacterota 
         ### Prevotellaceae family in order Bacteriodales

         
         
         # order the table so it matches the 
         # make otu column 
         # grab correct order 
         target<-tree.short$tip.label
         head(otu_log2)
         dim(otu_log2)
         otu_log2<-otu_log2[match(target, otu_log2$Family),]
         dim(otu_log2)
         # rm out column
         # makes sure no space in col names
         n <- otu_log2$Family
         otu_log2<-subset(otu_log2, select = c( -Family))
         row.names(otu_log2) <- n
         
         colnames(otu_log2)
         head(otu_log2)
         dim(otu_log2)
         tree.short
         #row.names(out_log2.1) <-f
       
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
         windows(10,10)
         heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))
         #dev.off()
         
         #make a plot without phylogeny 
         #make a dendrogram              
         row_dendro = as.dendrogram(hclust(dist((m))))
         
         ############  make the heatmap  plot #################33
         windows(10,10)
         heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
         
         
         
         

         
         
#____________________________________________#


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


##### Phyloseq and exploring  taxonomy ####

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


######------  import percent abundance into phyloseq for figure ----- #####

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




