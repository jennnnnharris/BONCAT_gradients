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
library(ape)
library(readxl)
#install.packages(ggtree)
#install.packages("TreeTools")
library('TreeTools')

library(NatParksPalettes)
library(MicEco)
library(ggtree)


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

library(NatParksPalettes)
names(NatParksPalettes)

natparks.pals(name="Arches",n=8,type="discrete")
natparks.pals(name="Arches2",n=4,type="discrete")

#get colors from national park pal

#DNA bulk soil
"#b46db3"
#DNA rhizo 
"#3a1f46"

#total rhizo
"#8fcafd"
#total root
"#4499f5"
#total nodule
"#0c62af"

#active rhizo
"#f0ac7d"
#active root
"#cd622e"
#active nodule
"#993203"

mycols6<-c("#8fcafd",  "#4499f5" ,  "#0c62af" ,   "#f0ac7d", "#cd622e" , "#993203")

mycols8 <- c("#b46db3", "#3a1f46", "#8fcafd",  "#4499f5" ,  "#0c62af" ,   "#f0ac7d", "#cd622e" , "#993203")
  

######-------------import Asvs data -----------------##############
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")


### Import Data ###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
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
#min.s<-min(rowSums(asvs.t))

### Rarefy to obtain even numbers of reads by sample ###
#set.seed(336)
#asvs.r<-rrarefy(asvs.t, min.s)
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
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

#test it worked
#sample_names(ps)
print(ps)
# 12855 taxa

# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that aren; in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
# 12855 taxa

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

#svg(file="figures/16s/diversity.svg",width = 10, height=4 )
#windows()
#plot_richness(ps, "Fraction", measures = c("Observed","Shannon", "Simpson", "InvSimpson")) +
#  geom_boxplot(aes(fill = "Fraction")) + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578"))
#dev.off()

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

rich$Fraction <-factor(rich$Fraction, levels = c("Total_DNA", "Total_Cells", "BONCAT_Active"))





# total + active number otus
svg(file="figures/16s/shannon_divserity.svg",width = 10, height=6 )
windows(width = 6, height=7)
rich %>%
  ggplot(aes(x=Compartment, y=Shannon,  col= Fraction))+
  geom_boxplot() +
  scale_color_manual(values=mycols8[c(4,7,1)]) +
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")+
  facet_wrap(~Fraction, scales = "free_x")+
  scale_x_discrete(drop = TRUE) +
  ylab("Shannon Diversity")+
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
#df<-rich %>%
#  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="SYBR") %>%
#  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
#  filter(Compartment=="Rhizosphere")
#m1<-lm(Observed ~ BONCAT,  data = df)
#summary(m1)
# sybr and DNA are different
#BONCATSYBR   -267.00      59.18  -4.512  0.00197 ** 


# bulk soil dna verse rhizosphere DNA
#df<-rich %>%
#  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="DNA") %>%
#  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
#  filter(Compartment=="Rhizosphere"| Compartment=="Bulk_Soil")
#m1<-lm(Observed ~ Compartment,  data = df)
#summary(m1)
## rhizosphere DNA verse bulk soil dna are not different
#CompartmentRhizosphere  -100.60      96.51  -1.042     0.32    



# nodule verse endosphere boncart
#df<-rich %>%
#  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="POS") %>%
#  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
#  filter(Compartment=="Nodule"| Compartment=="Roots")
#m1<-lm(Observed ~ Compartment,  data = df)
#summary(m1)

# nodule verse endosphere total

#df<-rich %>%
#  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
#  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
#  filter(Compartment=="Nodule"| Compartment=="Roots")
#m1<-lm(Observed ~ Compartment,  data = df)
#summary(m1)






#####DIVERSITY Stats #########
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
####### are taxa in the root present in the bulk soil ##########
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


################################# BARPLOT of active verse total abundance #####################
# plan
#phyla level
# phyla level %>% sum by phyla level %>% average replicates %>% rhizo 
                                                            #%>% root
                                                            #%>% nod
                                                            # view abundance total active

#otu level
# otu active  %>% average replicates %>% just rhizosphere %>% top 25 otus active %>%


#phyla level
# filter for total and active
sample_data(ps)
df<-subset_samples(ps, Compartment!="ctl"& Compartment!="Bulksoil" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df<-as.data.frame(t(otu_table(df)))
dim(df)
# 5771 taxa
df<-cbind(taxon,df)

# aggregate by phyl level 
phya<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29,
                      C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,  C1R.SYBR_S16,  C2E.SYBR_S22,  C2N.SYBR_S15,  C2R.SYBR_S17,
                      C5E.SYBR_S23, C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25,  C7R.SYBR_S19 ) ~ Phyla, data = df, FUN = sum, na.rm = TRUE)
dim(phya)
phya
# 44 ohyla

# average by rep
#nodule
phya$nodule.total.mean <-   rowMeans(phya %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
phya$nodule.total.log <- rowMeans(phya %>% dplyr::select(contains("N.SYBR"))) %>% log10() %>%  glimpse()
phya$nodule.xbcat.mean <-   rowMeans(phya %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
phya$nodule.xbcat.log <- rowMeans(phya %>% dplyr::select(contains("N.POS"))) %>% log10() %>%  glimpse()
# roots
phya$root.total.mean <-   rowMeans(phya %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
phya$root.total.log <- rowMeans(phya %>% dplyr::select(contains("E.SYBR"))) %>% log10() %>%  glimpse()
phya$root.xbcat.mean <-   rowMeans(phya %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
phya$root.xbcat.log <- rowMeans(phya %>%  dplyr::select(contains("E.POS"))) %>% log10() %>%  glimpse()

#rhizo
phya$rhizo.total.mean <-   rowMeans(phya %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
# add 1 to everything so there are no negative number
phya$rhizo.total.log<-(phya$rhizo.total.mean +1 ) %>% log10()

phya$rhizo.bcat.mean <-   rowMeans(phya %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
# add 1 to everything so there are no negative number
phya$rhizo.bcat.log<-(phya$rhizo.bcat.mean +1 ) %>% log10()

phya[phya== "-Inf"] <- 0
colnames(phya)
#plot to check
ggplot(phya, aes(x=Phyla , y=rhizo.total.log)) + 
  geom_bar(stat="identity", position="identity") +
  coord_flip()


# sep values, make group column, make total values negative
bcat<-phya%>% select(contains("bcat"))
bcat<-cbind(phya$Phyla, bcat)
bcat$group <- rep("ACTIVE", length(bcat$rhizo.bcat.mean))
dim(bcat)
bcat
total<-phya%>% select(contains("total"))
total<-cbind(phya$Phyla, total)
total$group <- rep("TOTAL", length(total$nodule.total.mean))
total<-total%>%mutate_if(is.numeric, funs(. *-1))


# rename colums
colnames(total)
colnames(bcat)
n<-c("Phyla", "nodule.mean", "nodule.log", "root.mean", "root.log", "rhizo.mean", "rhizo.log",  "group")
colnames(bcat)<-n
colnames(total)<-n
total
bcat

phyladf<-rbind(total, bcat)


#plot to check
windows(6,6)
ggplot(phyladf, aes(x=Phyla , y=rhizo.log, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("#ab8af2",
                             "#4c4b4d"))+
  xlab("phyla")+
  theme_bw()+
  coord_flip()

windows(6,6)
ggplot(phyladf, aes(x=Phyla , y=root.log, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("#f0c413",
                             "#4c4b4d"))+
  xlab("phyla")+
  theme_bw()+
  coord_flip()

windows(6,6)
ggplot(phyladf, aes(x=Phyla , y=nodule.log, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("#f75477",
                             "#4c4b4d"))+
  xlab("phyla")+
  theme_bw()+
  coord_flip()

# BAR CHART with TOP 50 OTUS
#nodule
df$nodule.total.mean <-   rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
df$nodule.total.log <- rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>% log10() %>%  glimpse()
df$nodule.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$nodule.xbcat.log <- rowMeans(df %>% dplyr::select(contains("N.POS"))) %>% log10() %>%  glimpse()
# roots
df$root.total.mean <-   rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
df$root.total.log <- rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% log10() %>%  glimpse()
df$root.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df$root.xbcat.log <- rowMeans(df %>%  dplyr::select(contains("E.POS"))) %>% log10() %>%  glimpse()

#rhizo
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
# add 1 to everything so there are no negative number
df$rhizo.total.log<-(df$rhizo.total.mean +1 ) %>% log10()

df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
# add 1 to everything so there are no negative number
df$rhizo.bcat.log<-(df$rhizo.bcat.mean +1 ) %>% log10()

df[df== "-Inf"] <- 0
colnames(df)

# select top 50 otus rhizosphere
otu<-df%>% select(-contains("POS")) %>% select(-contains("SYBR")) %>% select(-contains("nodule")) %>% select(-contains("root"))
otu<-otu[order(otu$rhizo.total.mean, decreasing = TRUE),]
rhizo<-otu[1:50,]

# sep values, make group column, make total values negative
bcat<-rhizo%>% select(contains("bcat"))
bcat<-cbind(rhizo$Genus, bcat)
bcat$group <- rep("ACTIVE", length(bcat$rhizo.bcat.mean))
dim(bcat)
bcat
total<-rhizo%>% select(contains("total"))
total<-cbind(rhizo$Genus, total)
total$group <- rep("TOTAL", length(total$`rhizo$Genus`))
total<-total%>%mutate_if(is.numeric, funs(. *-1))

# names
colnames(rhizo)
rhizo$otu <- rownames(rhizo)

#plot to check
ggplot(rhizo, aes(x=otu, y=rhizo.bcat.log)) + 
  geom_bar(stat="identity", position="identity") +
  coord_flip()

# rename colums
colnames(total)
colnames(bcat)
n<-c("Genus", "rhizo.mean", "rhizo.log",  "group")
colnames(bcat)<-n
colnames(total)<-n
total
bcat

rhizodf<-rbind(total, bcat)
rhizodf$otu <- rownames(rhizo)

#plot to check
ggplot(rhizodf, aes(x=otu , y=rhizo.log, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("#ab8af2",
                             "#4c4b4d"))+
  xlab("top 50 otus")+
  theme_bw()+
  coord_flip()
dim(df)
##### top roots
otu<-df%>% select(-contains("POS")) %>% select(-contains("SYBR")) %>% select(-contains("nodule")) %>% select(-contains("rhizo"))
otu<-otu[order(otu$root.total.mean, decreasing = TRUE),]
root<-otu[1:50,]
root

# sep values, make group column, make total values negative
bcat<-root%>% select(contains("bcat"))
bcat<-cbind(root$Genus, bcat)
bcat$group <- rep("ACTIVE", length(bcat$root.xbcat.mean))
dim(bcat)
bcat
total<-root%>% select(contains("total"))
total<-cbind(root$Genus, total)
total$group <- rep("TOTAL", length(total$`root$Genus`))
total<-total%>%mutate_if(is.numeric, funs(. *-1))
total
# names
colnames(rhizo)
root$otu <- rownames(root)

#plot to check
ggplot(root, aes(x=otu, y=root.xbcat.log)) + 
  geom_bar(stat="identity", position="identity") +
  coord_flip()

# rename colums
colnames(total)
colnames(bcat)
n<-c("Genus", "root.mean", "root.log",  "group")
colnames(bcat)<-n
colnames(total)<-n
total
bcat

rootdf<-rbind(total, bcat)
rootdf$otu <- rownames(root)

#rootdf[order(rootdf$root.mean),]


#plot to check
ggplot(rootdf, aes(x=otu , y=root.log, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("#f0c413",
                             "#4c4b4d"))+
  xlab("top 50 otus")+
  theme_bw()+
  coord_flip()

# nod
dim(df)
df
##### top nodule
otu<-df%>% select(-contains("POS")) %>% select(-contains("SYBR")) %>% select(-contains("root")) %>% select(-contains("rhizo"))
otu<-otu[order(otu$nodule.total.mean, decreasing = TRUE),]
nod<-otu[1:50,]
nod

# sep values, make group column, make total values negative
bcat<-root%>% select(contains("bcat"))
bcat<-cbind(root$Genus, bcat)
bcat$group <- rep("ACTIVE", length(bcat$root.xbcat.mean))
dim(bcat)
bcat
total<-root%>% select(contains("total"))
total<-cbind(root$Genus, total)
total$group <- rep("TOTAL", length(total$`root$Genus`))
total<-total%>%mutate_if(is.numeric, funs(. *-1))
total
# names
colnames(rhizo)
root$otu <- rownames(root)

#plot to check
ggplot(root, aes(x=otu, y=root.xbcat.log)) + 
  geom_bar(stat="identity", position="identity") +
  coord_flip()

# rename colums
colnames(total)
colnames(bcat)
n<-c("Genus", "nod.mean", "nod.log",  "group")
colnames(bcat)<-n
colnames(total)<-n
total
bcat

noddf<-rbind(total, bcat)
noddf$otu <- rownames(root)
noddf

#plot to check
ggplot(noddf, aes(x=otu , y=nod.log, fill=group)) + 
  geom_bar(stat="identity", position="identity") +
  scale_fill_manual(values=c("#f75477",
                             "#4c4b4d"))+
  xlab("top 50 otus")+
  theme_bw()+
  coord_flip()

######################### ANCOM ###############
#library(BiocManager)
#BiocManager::install("microbiome")

knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)


library(ANCOMBC)
library(tidyverse)
library(microbiome)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

##########example data

data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)

# subset to baseline
tse = tse[, tse$time == 0]

# Re-code the bmi group
tse$bmi = recode(tse$bmi_group,
                 obese = "obese",
                 severeobese = "obese",
                 morbidobese = "obese")
# Subset to lean, overweight, and obese subjects
tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]

# Note that by default, levels of a categorical variable in R are sorted 
# alphabetically. In this case, the reference level for `bmi` will be 
# `lean`. To manually change the reference level, for instance, setting `obese`
# as the reference level, use:
tse$bmi = factor(tse$bmi, levels = c("obese", "overweight", "lean"))
# You can verify the change by checking:
# levels(sample_data(tse)$bmi)

# Create the region variable
tse$region = recode(as.character(tse$nationality),
                    Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                    CentralEurope = "CE", EasternEurope = "EE",
                    .missing = "unknown")

# Discard "EE" as it contains only 1 subject
# Discard subjects with missing values of region
tse = tse[, ! tse$region %in% c("EE", "unknown")]

print(tse)


#########################VENN DIAGRAM are most taxa generalists or specialists#################################
### don't include taxa that are only in 1 samples and less than 50 reads total
sample_data(ps)
# at least in 2 samples min reads is 10
df<-subset_samples(ps, Fraction=="BONCAT_Active")
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df
# 852 taxa

df<-prune_taxa(taxa_sums(df) > 0, df)
# subset for active
#df<-prune_taxa(taxa_sums(df) > 5, df)
any(taxa_sums(df) == 0)

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
head(df)
habitat<-rep(NA, 2621)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>1]="rhizo specailist"
habitat[df$nodule==0 & df$roots>1 & df$rhizo==0]="root specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
#habitat
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
head(df)
nodule<-rownames(df[df$nodule==1,])
roots<-rownames(df[df$roots==1,])
rhizo<-rownames(df[df$rhizo==1,])
df
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
# otu level
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/venndiagram")
svg(file="OTU_level_active_venn.svg",width = 6, height=6 )


ggvenn(
  x, 
  fill_color = mycols6[c(4,5,6)],
  stroke_size = 2, set_name_size = 9, text_size = 6, digits = 1, fill_alpha=.8
  #auto_scale = TRUE
 ) +
  ggtitle("Active Otus")
dev.off()
####genus level??




###########-------- total cells 
df<-subset_samples(ps, Fraction=="Total_Cells")
# at least in 2 samples min reads is 10
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df
# 1342 taxa
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

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/venndiagram")
svg(file="OTU_level_total_venn.svg",width = 6, height=6 )
ggvenn(
  x, 
  fill_color = mycols6[c(1,2,3)],
  stroke_size = 2, set_name_size = 9, text_size = 6, digits = 1, fill_alpha=.8
  #auto_scale = TRUE
) +
  ggtitle("Total Viable Otus")
dev.off()

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
windows(6,6)
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
#asvs.bray<-vegdist(otu_table(ps.a), method = "bray")
#ps.a


# Perform PCoA analysis of BC distances #
#asvs.pcoa <- cmdscale(asvs.bray, k=(42-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
#asvs.p <- asvs.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
#asvs.eig<-asvs.pcoa$eig
#perc.exp<-asvs.eig/(sum(asvs.eig))*100
#pe1<-perc.exp[1]
#pe2<-perc.exp[2]

#plant verse soil
#metadat
# subset

#p<-otus.p[which(metadat$BONCAT != "ctl",)]

# subset metadata
#as.factor(metadat$Compartment)
#as.factor(metadat$Fraction)
#as.factor(metadat$compartment_BCAT)
#levels(as.factor(metadat$compartment_BCAT))
#levels(as.factor(metadat$Fraction))

# colors :)
#gold <- "#FFB000"
#purple <- "#785EF0"
#blue <- "#739AFF"
#pink <- "#DC267F"

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/pcoa/Pcoa_ASVS_plantvsoil_raw.svg",width = 7, height=6 )

#windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
#ordiplot(asvs.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
#points(asvs.p, col=c("black"),
    #   pch=c(22, 21,22,24, 23)[as.factor(metadat$Fraction)],
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

#####PHYLOGENIC SIGNAL STARTS HERE#########
####import feature table

######import asvs
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)
# filtered asvs tree to match with tree
asvs.filtered <- read.table("trees/SEPPoutput/feature-table.tsv", sep="\t", header=T, row.names = 1 )


## Transpose ASVS table ##
asvs.t <- t(asvs.filtered)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs.t<-asvs.t[order(row.names(asvs.t)),]

###--- recode metadata----- #
metadat<- metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule')) %>%
  mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl')) %>%
  select(-BONCAT) %>%
  mutate(compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction)) %>%
  glimpse()

### put into phyloseq
##------make phyloseq object with rarefied data -------#

asvs.phyloseq<- (asvs.t)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(asvs.filtered)
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
# 13766 taxa taxa



##------filter for active taxa -------#
glimpse(asvs.t)
dim(asvs.t)
ps1<-subset_samples(ps, Fraction=="BONCAT_Active")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 3319 taxa 
# at least in 2 samples min reads is 10
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 10)
ps1
any(taxa_sums(ps1) == 0)
# 850 taxa
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))

## summarise replicates
head(df)
#rhizo
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
#nodule
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
# roots
df$root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
colnames(df)
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
dim(df)

#### keep top 50 OTUS from each environment
head(df)
nod<-df[order(df$nodule, decreasing= TRUE),]
nod<-row.names(nod[1:50,])

head(df)
root<-df[order(df$root, decreasing= TRUE),]
root<-row.names(root[1:50,])

head(df)
rhizo<-df[order(df$rhizo, decreasing= TRUE),]
rhizo<-row.names(rhizo[1:50,])

otus<- unique(c(nod, root, rhizo))
length(otus)
# grab otus
df<-df[match(otus, df$otu),]
df<-df%>% filter(otu!="Others")
head(df)
dim(df)

df<-df[,c(4,3,2,1)]
# who are the taxa?
taxon$otu <- row.names(taxon)
lil_taxon<-left_join(df, taxon)
head(lil_taxon)
# make a list of genus
# if genus is missing replace with higher taxa
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Family[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Order[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Class[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Phyla[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Domain[which(lil_taxon$Genus == '')]

lil_taxon$Genus
##
lil_taxon$Genus[which(lil_taxon$Genus == ' g__')] <- lil_taxon$Family[which(lil_taxon$Genus == ' g__')]
lil_taxon$Genus[which(lil_taxon$Genus == ' f__')] <- lil_taxon$Order[which(lil_taxon$Genus == ' f__')]
lil_taxon$Genus[which(lil_taxon$Genus == ' o__')] <- lil_taxon$Class[which(lil_taxon$Genus == ' o__')]

## also if genus is just number replace with higher taxa
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Family[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Order[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Class[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Phyla[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Domain[which(lil_taxon$Genus == '')]

lil_taxon$Genus[grep("017", lil_taxon$Genus )] <- lil_taxon$Family[grep("017", lil_taxon$Genus) ]
lil_taxon$Genus[grep("21", lil_taxon$Genus )] <- lil_taxon$Family[grep("21", lil_taxon$Genus) ]
lil_taxon$Genus[grep("14", lil_taxon$Genus )] <- lil_taxon$Order[grep("14", lil_taxon$Genus) ]

lil_taxon$Genus<-paste0(1:length(lil_taxon$Genus), lil_taxon$Genus)

###############Import tree##########
#tree <- read.tree("2022.10.phylogeny.asv.nwk")
#tree <- read.tree("2022.10.phylogeny.id.nwk")
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree <- read.tree("trees/SEPPoutput/tree.nwk" )
tree$tip.label
taxon$Species
tree
# make tree smaller
asvs_remove<-setdiff(tree$tip.label, df$otu) #asvs we don't want
tree.short<-ape::drop.tip(tree, asvs_remove,  collapse.singles = TRUE,) # remove asvs we don't need
tree.short

# grab correct order 
target<-tree.short$tip.label
df<-df[match(target, df$otu),]
dim(df)
head(df)
# rm Genus column from df
df<-subset(df, select = c( -otu))
head(df)

############### make matrix
m<- df
m
m<-log10(m)
m[m== "-Inf"] <- 0
summary(m)
m[m==0] <- NA
m
#m$genus <- lil_taxon$Genus
m


#### make tip labels genus
tree.short$tip.label <- lil_taxon$Genus
tree.short

##### otu heatmap  ####  
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/heatmap")
#######heat map active top 150 otus
#m$genus <- lil_taxon$Genus
#m<-m[,c(4,1,2,3)]
#m<-m[,c(2,3,4)]
row.names(m) <- lil_taxon$Genus
m
svg(file="OTU_level_active_heatmap.svg",width = 8, height=8 )
windows(8,8)
p <- ggtree(tree.short, branch.length = .5) + 
  #xlim_tree(5) +
  geom_tiplab(size=3, align=TRUE, linesize=1, offset = 0) + 
  theme_tree2()
  #geom_cladelab()+
  #geom_nodelab( size=2, color="green", nudge_x = 1)+
 # nodelabels(node=tree.short$edge)
p
#p
gheatmap(p, m, 
         colnames=FALSE,
          legend_title="active taxa", offset = .5) +
  scale_x_ggtree() + 
  scale_fill_gradient(low="#fffaa2", high = "#bb0000", aesthetics = "fill", na.value = "white",
                      name="Abundance in Active")+
  ggtitle("Top 50 ASVS in each compartment")
 
######## otu level heatmap for total ##########
# filter for active taxa
sample_data(ps)
ps1<-subset_samples(ps, Fraction=="Total_Cells")
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 4187 taxa 
# at least in 2 samples min reads is 10
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 10)
ps1
any(taxa_sums(ps1) == 0)
# 1344 taxa
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))

## summarise replicates
head(df)
#rhizo
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB"))) %>%   glimpse()
#nodule
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB"))) %>%   glimpse()
# roots
df$root <-   rowMeans(df %>% dplyr::select(contains("E.SYB"))) %>%   glimpse()
colnames(df)
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
dim(df)
head(df)
# maybe keep taxa that are in the top 10% ??


#### keep top 50 OTUS from each environment
head(df)
nod<-df[order(df$nodule, decreasing= TRUE),]
nod[1:50,] # only 48 are positive 
nod<-row.names(nod[1:50,])

head(df)
root<-df[order(df$root, decreasing= TRUE),]
root[1:50,]
root<-row.names(root[1:50,])

head(df)
rhizo<-df[order(df$rhizo, decreasing= TRUE),]
rhizo<-row.names(rhizo[1:50,])

otus<- unique(c(nod, root, rhizo))
length(otus)
# grab otus
df<-df[match(otus, df$otu),]
df<-df%>% filter(otu!="Others")
head(df)
dim(df)

df<-df[,c(4,3,2,1)]

#
# who are the taxa?
taxon$otu <- row.names(taxon)
lil_taxon<-left_join(df, taxon)
head(lil_taxon)
# make a list of genus
# if genus is missing replace with higher taxa
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Family[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Order[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Class[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Phyla[which(lil_taxon$Genus == '')]
lil_taxon$Genus[which(lil_taxon$Genus == '')] <- lil_taxon$Domain[which(lil_taxon$Genus == '')]

lil_taxon$Genus
##
lil_taxon$Genus[which(lil_taxon$Genus == ' g__')] <- lil_taxon$Family[which(lil_taxon$Genus == ' g__')]
lil_taxon$Genus[which(lil_taxon$Genus == ' f__')] <- lil_taxon$Order[which(lil_taxon$Genus == ' f__')]
lil_taxon$Genus[which(lil_taxon$Genus == ' o__')] <- lil_taxon$Class[which(lil_taxon$Genus == ' o__')]

## also if genus is just number replace with higher taxa
lil_taxon$Genus[grep("017", lil_taxon$Genus )] <- lil_taxon$Family[grep("017", lil_taxon$Genus) ]
lil_taxon$Genus[grep("21", lil_taxon$Genus )] <- lil_taxon$Family[grep("21", lil_taxon$Genus) ]
lil_taxon$Genus[grep("14", lil_taxon$Genus )] <- lil_taxon$Order[grep("14", lil_taxon$Genus) ]
lil_taxon$Genus<-paste0(1:length(lil_taxon$Genus), lil_taxon$Genus)

## Import tree ##
#tree <- read.tree("2022.10.phylogeny.asv.nwk")
#tree <- read.tree("2022.10.phylogeny.id.nwk")
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree <- read.tree("trees/SEPPoutput/tree.nwk" )
tree$tip.label
taxon$Species

# make tree smaller
asvs_remove<-setdiff(tree$tip.label, df$otu) #asvs we don't want
tree.short<-ape::drop.tip(tree, asvs_remove,  collapse.singles = TRUE,) # remove asvs we don't need
tree.short

# grab correct order 
target<-tree.short$tip.label
df<-df[match(target, df$otu),]
dim(df)
head(df)

# rm Genus column from df
df<-subset(df, select = c( -otu))

# make tip names genus
tree.short$tip.label <- lil_taxon$Genus

###### make matrix
m<-df
#row.names(m)
#m<-structure(m)
# summary stats
# i think this matrix needs a log transform
m<-log10(m)
m[m== "-Inf"] <- 0
m
summary(m)
m[m==0] <- NA
#make row name shte genus
row.names(m) <- lil_taxon$Genus

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/heatmap")
#######heat map active top 150 otus
svg(file="OTU_level_total_heatmap.svg",width = 8, height=8 )

windows(8,8)
p <- ggtree(tree.short, branch.length = .1) + 
  #xlim_tree(5) +
  geom_tiplab(size=3, align=TRUE, linesize=1, offset = -.5) + 
  theme_tree2()

p
gheatmap(p, m, 
         colnames=FALSE, legend_title="active taxa") +
  scale_x_ggtree() + 
  # scale_color_gradient(l
  scale_fill_gradient(low= mycols6[1], high = "#0c328a", aesthetics = "fill", na.value = "white",
                      name="Abundance in TOTAL")+
  ggtitle("Top 50 ASVS in each compartment")
dev.off()




###### graveyard HEAT MAP group my habitat use type ######

#df
habitat<-rep(NA, 348)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="rhizospecailist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="rootspecialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodulespecialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizoandrootgeneralist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plantgeneralist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hypergeneralist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="rhizoandnodulegeneralist"

habitat
unique(habitat)
df$habitat <- habitat
df_habitat <-df
head(df_habitat)

# grab correct order 
target<-tree$tip.label
df_habitat<-df_habitat[match(target, df_habitat$Genus),]
glimpse(df_habitat)

#who are they?
n<-df_habitat%>% filter(habitat=="rootspecialist") %>% select(Genus) 
n
n<-df_habitat%>% filter(habitat=="rhizospecailist") %>% select(Genus) 
n<-df_habitat%>% filter(habitat=="rhizoandrootgeneralist") %>% select(Genus) 
n<-df_habitat%>% filter(habitat=="plantgeneralist") %>% select(Genus) 
n<-df_habitat%>% filter(habitat=="hypergeneralist") %>% select(Genus) 

n
unique(df_habitat$habitat)
############### make matrix
df<-df_habitat %>% select(rhizo,roots, nodule)
n<-row.names(df)
m<-as.matrix(df)
m<-structure(m)


## i think this matrix needs a log transform
m<-log10(m)
head(m)
m[m== "-Inf"] <- 0
head(m)
summary(m)
# min = 0
# max = 5.099
glimpse(m)
tree

#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/heatmap2.svg", width = 10, height=10 )
# And make the plot with phylogeny
windows(16,12)
heatmap.phylo(x = m, Rowp = tree, Colp = as.phylo(as.hclust(col_dendro)))+
  
  dev.off()
# change order
#tr<-read.tree(text =  "((Nodule1:3,Nodule2:3, Nodule3:3,Nodule4:3):2, (Root1:3,Root2:3,Root3:3,Root4:3):2, (Rhizosphere1:3, Rhizosphere2:3,Rhizosphere3:3, Rhizosphere4:3):2);")
heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)
dev.off
# no phylogeny
#make a dendrogram              
row_dendro = as.dendrogram(hclust(dist((m))))
windows(10,10)
heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
#legend(legend(x="right", legend=c("min", "med", "max"),fill=

dev.off


















############################################################### phylogeny graveyard
# add taxa names
#n<-row.names(otu)
#otu<-mutate(otu, OTU=n)
#n<-row.names(taxon)
#taxon<- mutate(taxon, OTU=n)
#otu<-left_join(otu, taxon)


#otu
# put g__ in front
#otu$Genus <- gsub('\\ ' , '', otu$Genus)
#otu$Genus<-paste0( "g__",otu$Genus)
#otu$Family <- gsub('\\ ' , '', otu$Family)
##otu$Family<-paste0( "f__",otu$Family)
#otu$Order <- gsub('\\ ' , '', otu$Order)
#otu$Order<-paste0( "o__",otu$Order)
#otu$Class <- gsub('\\ ' , '', otu$Class)
#otu$Class<-paste0( "c__",otu$Class)
#otu$Phyla <- gsub('\\ ' , '', otu$Phyla)
#otu$Phyla<-paste0( "p__",otu$Phyla)


# if genus is missing put family name, if still missing put order, class, domain, etc. 
#otu$Genus[which(otu$Genus == 'g__')] <- otu$Family[which(otu$Genus == 'g__')]
#otu$Genus[which(otu$Genus == 'f__')] <- otu$Order[which(otu$Genus == 'f__')]
#otu$Genus[which(otu$Genus == 'o__')] <- otu$Class[which(otu$Genus == 'o__')]
#otu$Genus[which(otu$Genus == 'c__')] <- otu$Phyla[which(otu$Genus == 'c__')]
#otu$Genus[which(otu$Genus == 'p__')] <- otu$Domain[which(otu$Genus == 'p__')]

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
#length(unique(otu$Genus))
# 499 bacteria unique taxa groups


###changing my data's names so they match the tree ####################

# put order as same for type III 
#otu$Genus[grepl("g__type_III", otu$Genus )] <- otu$Order[grepl("g__type_III", otu$Genus) ]
#g__Tepidisphaeraceae" is a family
#otu$Genus[grepl("g__Tepidisphaeraceae", otu$Genus)] <- otu$Family[grepl("g__Tepidisphaeraceae", otu$Genus)]
# switch "g__Nitrospira_ to "g__Nitrospira_A" 
#otu$Genus <- sub("g__Nitrospira" ,  "g__Nitrospira_A"   ,  otu$Genus )
# "g__Vicinamibacteraceae" is a family
#otu$Genus[grepl("g__Vicinamibacteraceae", otu$Genus)] <-  otu$Family[grepl("g__Vicinamibacteraceae", otu$Genus )] 
#  "g__Vampirovibrionaceae"  is a family
#otu$Genus[grepl("g__Vampirovibrionaceae" , otu$Genus)] <-  otu$Family[grepl("g__Vampirovibrionaceae" , otu$Genus )] 
# "g__Vermiphilaceae"     is family
#otu$Genus[grepl( "g__Vermiphilaceae"   , otu$Genus)] <-  otu$Family[grepl(  "g__Vermiphilaceae" , otu$Genus )] 
# "g__Vampirovibrionales"  is a order
#otu$Genus[grepl( "g__Vampirovibrionales" , otu$Genus)] <-  otu$Order[grepl("g__Vampirovibrionales" , otu$Genus )] 
#### change g__Peribacteria" to  "g__Peribacter"
#otu$Genus <- sub("g__Peribacteria", "g__Peribacter" ,  otu$Genus )
## Parcubacteria is a phyla so make it a phyla
#otu$Genus[grepl("g__Parcubacteria", otu$Genus)] <-  otu$Phyla[grepl("g__Parcubacteria", otu$Genus )] 

### _Pedosphaeraceae is a family, call it a family
#otu$Genus[grepl("g__Pedosphaeraceae", otu$Genus)] <-  otu$Family[grepl("g__Pedosphaeraceae", otu$Genus )] 


### replace unknown family with order
#otu$Genus[grepl("g__Unknown_Family", otu$Genus)] <-  otu$Order[grepl("g__Unknown_Family", otu$Genus )] 

### Berkelbacteria is actual a phyla
#otu$Genus[grepl("g__Berkelbacteria", otu$Genus)] <-  otu$Phyla[grepl("g__Berkelbacteria", otu$Genus )] 

### Gracilibacteria is actual a phyla
#otu$Genus[grepl("g__Gracilibacteria", otu$Genus)] <-  otu$Phyla[grepl("g__Gracilibacteria", otu$Genus )] 

#"g__Chthoniobacteraceae" is an family, but only 1 genus in it, so I'm gonna call it that 
#otu$Genus <- sub("g__Chthoniobacteraceae", "g__Chthonomonas" ,  otu$Genus , ignore.case = TRUE)

#"g__Chthoniobacteraceae" is an order, but only 1 genus in it, so I'm gonna call it that 
#otu$Genus <- sub("g__Chthonomonadales", "g__Chthonomonas" ,  otu$Genus , ignore.case = TRUE)

#"g__#Entotheonaella" is an family, but only 1 genus in it, so I'm gonna call it that 
#otu$Genus <- sub("g__Entotheonellaceae", "g__Entotheonella" ,  otu$Genus )

#Fimbriimonadaceae is a family  , but only 1 genus in it, so I'm gonna call it that 
#otu$Genus <- sub("g__Fimbriimonadaceae", "g__Fimbriimonas" ,  otu$Genus )

# Fimbriimonadales is an order , but only 1 genus in it, so I'm gonna call it that 
#otu$Genus <- sub("g__Fimbriimonadales", "g__Fimbriimonas" ,  otu$Genus )

# Hyphomicrobium call it g__Hyphomicrobium_A
#otu$Genus <- sub("g__Hyphomicrobium", "g__Hyphomicrobium_A" ,  otu$Genus )

### change "g__Latescibacteraceae" to lascibacter
#otu$Genus <- sub("g__Latescibacterota", "g__Latescibacter" ,  otu$Genus )

### change "g__Latescibacterota" to lascibacter
#otu$Genus <- sub("g__Latescibacteraceae", "g__Latescibacter" ,  otu$Genus )

## change  to g__Massilia to "g__Massilibacterium"  because they are the same
#otu$Genus <- sub("g__Massilia", "g__Massilibacterium" ,  otu$Genus )

# Methyloligellaceae is a family so call it a family
#otu$Genus[grepl("Methyloligellaceae", otu$Genus)] <-  otu$Family[grep("Methyloligellaceae", otu$Genus )] 

# Microgenomatia is a phyla so call it a phyla
#otu$Genus[grepl("g__Microgenomatia", otu$Genus)] <-  otu$Phyla[grep("g__Microgenomatia", otu$Genus )] 

# Rokubacteriales is an order. call it an order
#otu$Genus[grepl("g__Rokubacteriales", otu$Genus)] <-  otu$Order[grep("g__Rokubacteriales", otu$Genus )] 


# replace linage with phyla
#otu$Genus[grepl("g__Lineage", otu$Genus)] <-  otu$Phyla[grep("g__Lineage", otu$Genus )] 

## anything with number replace with a higher taxonomic level
#otu$Genus[grep("[7]", otu$Genus)] <-  otu$Family[grep("[7]", otu$Genus )] 
#otu$Genus[grep("[7]", otu$Genus)] <-  otu$Order[grep("[7]", otu$Genus )] 
#otu$Genus[grep("[7]", otu$Genus)] <-  otu$Class[grep("[7]", otu$Genus )] 
#otu$Genus[grep("[7]", otu$Genus)] <-  otu$Phyla[grep("[7]", otu$Genus )] 

#otu$Genus[grep("[-]", otu$Genus)] <-  otu$Family[grep("[-]", otu$Genus )] 
#otu$Genus[grep("[-]", otu$Genus)] <-  otu$Order[grep("[-]", otu$Genus )] 
#otu$Genus[grep("[-]", otu$Genus)] <-  otu$Class[grep("[-]", otu$Genus )] 
#otu$Genus[grep("[-]", otu$Genus)] <-  otu$Phyla[grep("[-]", otu$Genus )] 

#otu$Genus[grep("[1]", otu$Genus)] <-  otu$Family[grep("[1]", otu$Genus )] 
#otu$Genus[grep("[1]", otu$Genus)] <-  otu$Order[grep("[1]", otu$Genus )] # A0
#otu$Genus[grep("[1]", otu$Genus)] <-  otu$Class[grep("[1]", otu$Genus )] 
#otu$Genus[grep("[1]", otu$Genus)] <-  otu$Phyla[grep("[1]", otu$Genus )] 
#otu$Genus[grep("[1]", otu$Genus)] <-  otu$Domain[grep("[1]", otu$Genus )] 

#otu$Genus[grep("[2]", otu$Genus)] <-  otu$Family[grep("[2]", otu$Genus )] 
#otu$Genus[grep("[2]", otu$Genus)] <-  otu$Order[grep("[2]", otu$Genus )] 
#otu$Genus[grep("[2]", otu$Genus)] <-  otu$Class[grep("[2]", otu$Genus )] 
#otu$Genus[grep("[2]", otu$Genus)] <-  otu$Phyla[grep("[2]", otu$Genus )] 
#otu$Genus[grep("[2]", otu$Genus)] <-  otu$Domain[grep("[2]", otu$Genus )] 

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
#otu$Genus <- sub("\\o__Absconditabacteriales_\\(SR1\\)", "o__Absconditabacteriales",  otu$Genus )

# swap "g__Clostridium_sensu_stricto_13"for "g__Clostridium" 
#otu$Genus <- sub("g__Clostridium_sensu_stricto_13", "g__Clostridium" ,  otu$Genus , ignore.case = TRUE)


# candidatus, because that just means it hasn't been cultured
#otu$Genus<- gsub('g__Candidatus_', 'g__', otu$Genus )

#change Armatimonadota to  "g__Armatimonas
#otu$Genus<- gsub('g__Armatimonadales', 'g__Armatimonas', otu$Genus )

## change g__Altererythrobacter to g__Altererythrobacter_D" 
#otu$Genus <- gsub('g__Altererythrobacter', 'g__Altererythrobacter_D', otu$Genus )

# for the groups with genusus that are just number use higher taxonomic classification
#otu$Genus[grepl('g__0', otu$Genus )] <- otu$Class[grepl('g__0', otu$Genus)]

# for the groups with uncultured as genus  use higher taxonomic classification
#otu$Genus[grepl('g__un', otu$Genus )] <- otu$Family[grepl('g__un', otu$Genus)]
#otu$Genus[grepl('f__un', otu$Genus )] <- otu$Order[grepl('f__un', otu$Genus)]
#otu$Genus[grepl('o__un', otu$Genus )] <- otu$Class[grepl('o__un', otu$Genus)]

# for the groups with subgroup as genus  use higher taxonomic classification
#otu$Genus[grepl('Subgroup', otu$Genus )] <- otu$Family[grepl('Subgroup', otu$Genus)]
#otu$Genus[grepl('f__Sub', otu$Genus )] <- otu$Order[grepl('f__Sub', otu$Genus)]
#otu$Genus[grepl('o__Sub', otu$Genus )] <- otu$Class[grepl('o__Sub', otu$Genus)]

# rename rhizobium group as such
#otu$Genus <- gsub('g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium' , 'g__Rhizobium', otu$Genus)

# rename burkholderia group as such
#otu$Genus <- gsub( 'g__Burkholderia-Caballeronia-Paraburkholderia' , 'g__Burkholderia', otu$Genus)

# replace Acidibacter with the order
#otu$Genus <- gsub('g__Acidibacter',  'o__Gammaproteobacteria_Incertae_Sedis' , otu$Genus)

# replace "g__Absconditabacteriales_(SR1)") with the order level
#otu$Genus[grepl('g__Abscond', otu$Genus )] <- otu$Order[grepl('Abscond', otu$Genus)]

##aggregate by genus

# how many genuses are there now?
length(unique(otu$Genus))
# 499 bacteria unique taxa groups

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

#------> skip this part if you want to group without phylogeny
 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
#bacteria
#tree = read.tree("bac120.tree")
#tree = read.tree("genus.nwk")
#tax = read.table("bac120_taxonomy.tsv")
#map = read.table("family.map")
# Set the working directory; ###
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/asv_level_output/greengenes")
### Import Data ###
#tree <- read.tree("2022.10.phylogeny.asv.nwk")
#tree <- read.tree("2022.10.phylogeny.id.nwk")
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree <- read.tree("trees/SEPPoutput/tree.nwk" )
tree$tip.label
taxon$Species

install.packages("rjson")
library("rjson")
placements<-fromJSON(file="trees/SEPPoutput/placements.json")
df<-as.data.frame(placements)


###making tip names genus  ## 
# put a g__ in front of all of these
#tree$tip.label<-paste0("g__", tree$tip.label)
#tree$tip.label
# remove random little /'-
#tree$tip.label<-gsub("\\'", "", tree$tip.label)
#length(tree$tip.label)
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

plot(tree.label, no.margin=TRUE,  cex = .5, show.node.label = TRUE)
nodelabels()
dev.off()
length(tree.label$tip.label)

tree.label$node.label


###tips to add###

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



##################################### heatmap genus level?? by phylogeny  ###############

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

## grab taxa that are annotated by greengenes 2 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
### Import Data ###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)
## Transpose ASVS table ##
asvs.t <- t(asvs.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs.t<-asvs.t[order(row.names(asvs.t)),]

###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

##---make phyloseq object with rarefied data -------#
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
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)
print(ps)
# 12855 taxa


##select active taxa 
# at least in 2 samples min reads is 50
ps1<-subset_samples(ps, Fraction=="BONCAT_Active")
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 50)
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))
dim(df)
# 726 taxa
####### aggregate to family level again
###### just doing family level again to keep it simple
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)
head(df)
#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#summarise
#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"

head(df)
#aggregate
df<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family

dim(df) # 140 families
#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS")))
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS")))
df$root <-   rowMeans(df %>% dplyr::select(contains("E.POS")))
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)
# for now rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>%   glimpse()
summary(df)
head(df)
df<-mutate(df, Family=otu) %>% select(., -otu)

dim(df)






habitat<-rep(NA, 138)
habitat[df$nodule==0 & df$root==0 & df$rhizo>0]="rhizospecailist"
habitat[df$nodule==0 & df$root>0 & df$rhizo==0]="rootspecialist"
habitat[df$nodule>0 & df$root==0 & df$rhizo==0]="nodulespecialist"
habitat[df$nodule==0 & df$root>0 & df$rhizo>0]="rhizoandrootgeneralist"
habitat[df$nodule>0 & df$root>0 & df$rhizo==0]="plantgeneralist"
habitat[df$nodule>0 & df$root>0 & df$rhizo>0]="hypergeneralist"
habitat[df$nodule>0 & df$root==0 & df$rhizo>0]="rhizoandnodulegeneralist"

habitat
unique(habitat)
df$habitat <- habitat
df


#who are they?
n<-df%>% filter(habitat=="rootspecialist") %>% select(Family) 
n
n<-df%>% filter(habitat=="rhizospecailist") %>% select(Family) 
n
n<-df%>% filter(habitat=="rhizoandrootgeneralist") %>% select(Family) 
n
n<-df%>% filter(habitat=="plantgeneralist") %>% select(Family) 
n
n<-df%>% filter(habitat=="hypergeneralist") %>% select(Family) 
n
n
unique(df_habitat$habitat)
############### make matrix
m<-df
m<-as.matrix(m[,c(3,2,1)])
dim(m)
#matrix needs a log transform
m<-m+1
m<-log10(m)
#m[m== "-Inf"] <- 0
#m[m== 0] <-NA
head(m)
summary(m)

heatmap(m)

ggplot(df, aes(X, Y, fill= Z)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdPu") +
  theme_ipsum()

# heatmap
# The mtcars dataset:
data <- as.matrix(mtcars)
data
# Default Heatmap
heatmap(data)
m

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

dat$rhizo.log <- rowMeans(otu %>% 
            dplyr::select(contains("Rhizo"))) %>% 
            log10() %>%
            glimpse() 

dat$root.mean <- rowMeans(otu %>% 
      dplyr::select(contains("Root"))) %>% 
      glimpse()

dat$root.log <- rowMeans(otu %>% 
                dplyr::select(contains("Root"))) %>% 
                log10() %>%
                glimpse()

dat$nodule.mean <- rowMeans(otu %>% 
      dplyr::select(contains("Nod"))) %>% 
      glimpse()

dat$nodule.log <- rowMeans(otu %>% 
                    dplyr::select(contains("Nod"))) %>% 
                    log10() %>%
                    glimpse()


dat$random <- rnorm(153, sd = 10)
dat$bm <- rTraitCont(tree.short)
#turn it into a data frame now?
dat <- as.data.frame(dat)
 
dat[dat== "-Inf"] <- 0 # change -inf into 0
head(dat)

#dim(tree.short)
tree.short
dat<-dat%>% select(contains("log"))

#We can combine phylogeny and traits into a phylo4d object.
p4d <- phylo4d(tree.short, dat)
windows(14,10)                               
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



#####--------heatmap  family level activity-------########
## grab taxa that are annotated by greengenes 2 
## because we are gonna use that tree
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
### Import Data ###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)
## Transpose ASVS table ##
asvs.t <- t(asvs.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs.t<-asvs.t[order(row.names(asvs.t)),]

###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

##---make phyloseq object with rarefied data -------#
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
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)
print(ps)
# 12855 taxa


##select active taxa 
# at least in 2 samples min reads is 50
ps1<-subset_samples(ps, Fraction=="BONCAT_Active")
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 50)
# 852 taxa
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))

# 726 taxa
####### agregate to the family level
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)
head(df)
#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"

head(df)
#aggregate
df<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family

dim(df) # 144 families
#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)
# for now rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>%   glimpse()
summary(df)
head(df)
df<-mutate(df, Family=otu) %>% select(., -otu)

###Import tree#
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/asv_level_output/greengenes")
#tree <- read.tree("2022.10.phylogeny.asv.nwk")
#tree <- read.tree("2022.10.phylogeny.id.nwk")

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")
tree
# 1057 tips
length(tree$tip.label) # look at the tip labels 
# modify tip labels
tree$tip.label <- paste0(" f__", tree$tip.label)

length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
# 69 intersect
mynames<-(setdiff(unique(df$Family), tree$tip.label))
length(mynames)

#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, df$Family) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE,  cex = .5)
# grab correct order 
target<-tree.short$tip.label
df<-df[match(target, df$Family),]
head(df)
# rm family column from df
df<-subset(df, select = c( -Family))
# place holder
m <- df
summary(m)
# i think this matrix needs a log transform
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]

# heat map
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/heatmap")
svg(file="Family_level_active_heatmap.svg",width = 8, height=8 )

windows(8,8)
p <- ggtree(tree.short, branch.length = .001) + 
  geom_tiplab(size=3, align=TRUE, linesize=0, offset = -.2) + 
  theme_tree2()+
  xlim_tree(1) 

#geom_cladelab()+
#geom_nodelab( size=2, color="green", nudge_x = 1)+
# nodelabels(node=tree.short$edge)
p
#p
gheatmap(p, m, 
         colnames=FALSE,
         legend_title="active taxa", offset = .5) +
  scale_x_ggtree() + 
  scale_fill_gradient(low="#fffaa2", high = "#bb0000", aesthetics = "fill", na.value = "white",
                      name="Abundance in Active")+
  ggtitle("Top 50 ASVS in each compartment")
dev.off()
#####--------heatmap family level total-------########
##select active taxa 
# at least in 2 samples min reads is 50
sample_data(ps)
ps1<-subset_samples(ps, Fraction=="Total_Cells")
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 50)
any(taxa_sums(ps1) == 0)
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))
dim(df)
# 971 taxa

####add family column
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)
head(df)
#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"
df$Family[grepl(" f__Bacillaceae_H", df$Family)] <- " Bacillaceae_H"
df$Family[grepl(" f__Streptomycetaceae", df$Family)] <-  " f__Streptomycetaceae"

#aggregate
df<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,  C1R.SYBR_S16,  C2E.SYBR_S22,
                    C2N.SYBR_S15,  C2R.SYBR_S17 ,  C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25,  C7R.SYBR_S19 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 145 families
#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.SYB"))) %>%   glimpse()

df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
#rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>%   glimpse()
df<-mutate(df, Family=otu) %>% select(., -otu)
head(df)

# grab phylogeny and Read in the tree file
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")
tree
# modify tip labels
tree$tip.label <- paste0(" f__", tree$tip.label)
length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
# 67 intersect
mynames<-(setdiff(unique(df$Family), tree$tip.label))
length(mynames)
mynames

#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, df$Family) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE,  cex = .5)
# grab correct order 
target<-tree.short$tip.label
df<-df[match(target, df$Family),]
head(df)
# rm family column from df
df<-subset(df, select = c( -Family))
# place holder
m <- df

summary(m)
# i think this matrix needs a log transform
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]

# heat map
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/heatmap")
svg(file="Family_level_total_heatmap.svg",width = 8, height=8 )

windows(8,8)
p <- ggtree(tree.short, branch.length = .001) + 
  geom_tiplab(size=3, align=TRUE, linesize=0, offset = -.2) + 
  theme_tree2()+
  xlim_tree(1) 
p

gheatmap(p, m, 
         colnames=FALSE,
         legend_title="active taxa", offset = .5) +
  scale_x_ggtree() + 
  scale_fill_gradient(low= "#8fcafd", high = "#0c328a", aesthetics = "fill", na.value = "white",
                      name="Abundance in TOTAL")+
    ggtitle("Families in Total Community")
dev.off()

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




