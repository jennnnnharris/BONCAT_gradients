# 16S analysis
# Jennifer Harris
# Last updated: May 7 2024

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

################## Load required libraries ############

#basic
library(tidyverse)
library(vegan)
library(reshape2)
library(scales)
library(data.table)


#Phyloseq and mbiome
library(phyloseq)
library(microbiome)
library(MicEco)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

#colors and patterns
library(RColorBrewer)

#tree
library(Heatplus)
library(ade4)
library(ape)
library('TreeTools')
library(ggtree)

# heatmaps 
library(gplots)

# venn diagrams
library(ggvenn)
library(ggplot2)
library(dplyr)
library(grid)

#Ancom
library(ANCOMBC)

## if trouble loading phyloseq, see: http://joey711.github.io/phyloseq/install
#source("http://bioconductor.org/biocLite.R")
#biocLite("Heatplus")
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("Heatplus")

#####Import data#####
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/")

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
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Viable_Cell', 'POS'='Active_Cell', 'ctl'= 'ctl'))
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

##------make phyloseq object with percent data -------#
asvs.phyloseq<- (asvs.t)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

# import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs.phyloseq), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps
# 12855 taxa 


#####remove plant contamination  ########
# NOTE: 
# we assigned taxonomy with the greengenes2 database
# greengenes2 doesn't have great taxonomy of chloroplast and mitchondria. 
# We selected taxa that were unassigned to any phyla in the the root or nodule.
# the taxonomy of these ASVS was determined with SILVA. 
# we used BLAST to identify any ASVS that were not idenitfied with SILVA
# we removed 48 ASVS that where mitochondria or chloroplasts. 

# Select unassigned Asvs that the only in the the roots and nodules
df<-subset_taxa(ps, Phyla=="" | Phyla == " p__")
df<-as.data.frame(t(otu_table(df)))
df<-select(df, contains ("E"))
df1<-as.data.frame(rowSums(df))
colnames(df1) <- "sum"
df1<-filter(df1, sum!=0)

df<-subset_samples(ps, Compartment=="Roots" | Compartment=="Nodule" )
df<-prune_taxa(taxa_sums(df) > 0, df)
remove<-subset_taxa(df, Domain=="Unassigned" |  Phyla=="" | Phyla==" p__"  ) 
remove

### These 5 ASVS were not mitochondria or chloroplasts, the remain 48 were
kp<-c("a49a51f3a3e3ea140206b10c5665cc13", "55c30bcbeacfedffa7aeb332600548b2" , "cfeae1df224b7e426ea125ab2bb824fc", "b260024f11a7d77d4f03e5ca2e239860", "f5211207035c7ea5b6f2ecfdad3765e1")
badtaxa<-taxa_names(remove)
badtaxa <- badtaxa[!(badtaxa %in% kp)]
length(badtaxa)

########## remove these ASVs 
alltaxa<-taxa_names(ps)
mytaxa <- alltaxa[!(alltaxa %in% badtaxa)]
ps<-prune_taxa(mytaxa, ps )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 12807 taxa
#### quick checks
# how many reads per samples after QC
n<-rowSums(otu_table(ps))
# number of ASVS per sample after QC
n<-otu_table(ps)
n<-rowSums(ifelse(n[]>0,1,0))

####PERCENT abundance figure: SUPPLEMENT####
# grab data
taxon<- as.data.frame(tax_table(ps))
df<-as.data.frame(otu_table(ps))

# make it percent
df<-(df/rowSums(df))*100
df<-as.data.frame(t(df))
df<-cbind(df, taxon)

# rename columns 
colnames(df)
length(n)
length(colnames(df))
n<-c("BEADS" ,    "Bulksoil.DNA.1"  , "Endo.BONCAT.1" , "Nodule.BONCAT.1" , "Nodule.Totalcells.1" ,"Rhizo.DNA.1",  "Rhizo.BONCAT.1" , "Rhizo.Totalcells.1" ,"Endo.BONCAT.2" ,  "Endo.Totalcells.2", 
     "Nodule.BONCAT.2" ,  "Nodule.Totalcells.2"  ,"Rhizo.DNA.2" ,   "Rhizo.BONCAT.2",   "Rhizo.Totalcells.2",  "Bulksoil.DNA.3",    "Endo.BONCAT.3",   "Endo.Totalcells.3" , "Nodule.BONCAT.3",   "Nodule.Totalcells.3" ,
     "Rhizo.DNA.3" ,  "Rhizo.BONCAT.3"  , "Rhizo.Totalcells.3" , "Bulksoil.DNA.4" ,    "Endo.BONCAT.4",    "Endo.Totalcells.4",  "Nodule.BONCAT.4" ,  "Rhizo.DNA.4",   "Rhizo.BONCAT.4",   "Rhizo.Totalcells.4" ,
     "BulkSoil.DNA.5",   "Endo.BONCAT.5",   "Endo.Totalcells.5",  "Nodule.BONCAT.5" ,  "Nodule.Totalcells.5" , "Rhizo.DNA.5",   "Rhizo.Totalcells.5" , "CTL"  ,     "Bulksoil.DNA.6"  ,  "Bulksoil.DNA.7"  ,  
     "Bulksoil.DNA.8",     "Bulksoil.DNA.9"  ,   "Domain"  ,      "Phyla"       ,  "Class"  ,       "Order"    ,     "Family"  ,      "Genus"    ,     "Species"   )
df1<-df
colnames(df1)<-n

# make rownames null
# summarize by phyla
df1<-aggregate(cbind(BEADS , Bulksoil.DNA.1  , Endo.BONCAT.1 , Nodule.BONCAT.1 , Nodule.Totalcells.1 , Rhizo.DNA.1,  Rhizo.BONCAT.1 , Rhizo.Totalcells.1 ,Endo.BONCAT.2 ,  Endo.Totalcells.2, 
                     Nodule.BONCAT.2 ,  Nodule.Totalcells.2  , Rhizo.DNA.2,   Rhizo.BONCAT.2,   Rhizo.Totalcells.2,  Bulksoil.DNA.3,    Endo.BONCAT.3,   Endo.Totalcells.3 , Nodule.BONCAT.3,   Nodule.Totalcells.3 ,
                     Rhizo.DNA.3 ,  Rhizo.BONCAT.3  , Rhizo.Totalcells.3 , Bulksoil.DNA.4 ,    Endo.BONCAT.4,    Endo.Totalcells.4,  Nodule.BONCAT.4 ,  Rhizo.DNA.4,   Rhizo.BONCAT.4,   Rhizo.Totalcells.4 ,
                     BulkSoil.DNA.5,   Endo.BONCAT.5,   Endo.Totalcells.5,  Nodule.BONCAT.5 ,  Nodule.Totalcells.5 , Rhizo.DNA.5,   Rhizo.Totalcells.5 , CTL  ,     Bulksoil.DNA.6  ,  Bulksoil.DNA.7  ,  
                     Bulksoil.DNA.8,     Bulksoil.DNA.9) ~ Phyla, data = df1, FUN = sum, na.rm = TRUE)

head(df1)
# summ row 1 and 2 b\c they are both unassigned taxa

row1<-df1[1,2:43]+ df1[2,2:43] 
# call empty phyla unassigned
row1<-c("Unassigned", row1)
# put in df
row1<-as.vector(row1)
df1[1,] <- row1
df1<-df1[c(1,3:53),]

# gather by sample
df1<-gather(df1, "sample", value, 2:43 )
head(df1)
#remove zeros
df1<-df1[df1$value!=0,]
head(df1)


# make really low abundance taxa other
df1$Phyla[df1$value<1] <- "other"
df1<-aggregate(cbind(value) ~ sample+Phyla, data = df1, FUN = sum, na.rm =TRUE)
head(df1)
df1<-df1[order(df1$sample),]
head(df1)


mycols18<- c( "#1F78B4","#A6CEE3","#E31A1C",  "#FB9A99", "#33A02C","#B2DF8A",  "#FF7F00",  "#FDBF6F", "#6A3D9A" , "#CAB2D6",
               "#B15928", "#FFFF99",  "#eb05db","#edceeb","#1a635a","#9ad6ce", "#969696", "#232423")

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/supplement")
svg(file="S3_barplot.svg",width = 12, height=10)
windows(12,10)
df1%>% 
  ggplot(aes(fill=Phyla, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values=mycols18) +
  #scale_fill_viridis(discrete = TRUE) +
  #ggtitle("Top phyla") +
  theme_bw(base_size = 12)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

dev.off()

####DIVERSITY  FIG 3####
# diversity is calculated on raw reads b/c many diversity metric use singleton to calculate diversity. 
# Normalizing data first can create an inaccurate estimate of diversity. 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

# Data wrangling fo rdiversity of active microbes in each fraction
rich<-cbind(rich, metadat)
rich<-as.data.frame(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))
rich %>%
arrange( -Observed)

rich<- rich %>%  filter(Plant!="NOPLANT", Fraction!="ctl")
rich$compartment_BCAT <-factor(rich$compartment_BCAT, levels = c("Bulk_SoilTotal_DNA", "RhizosphereTotal_DNA", "RhizosphereViable_Cell", "RhizosphereActive_Cell",
                                                                 "RootsViable_Cell"   ,  "RootsActive_Cell" , "NoduleViable_Cell" ,  "NoduleActive_Cell" ))          
rich$Fraction <-factor(rich$Fraction, levels = c("Total_DNA", "Viable_Cell", "Active_Cell"))
rich %>% group_by(rich$compartment_BCAT) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(rich$Compartment) %>% summarise(mean(Observed), sd(Observed))

mycols3 <- c("#bcd3e8",  "#282c55", "#fc8449")
#shannon
p1<-rich %>%
  ggplot(aes(x=Compartment, y=Shannon,  col= Fraction, fill=Fraction))+
  geom_boxplot(alpha=.5) +
  scale_color_manual(values=mycols3) +
  scale_fill_manual(values = mycols3)+
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")+
  facet_wrap(~Fraction, scales = "free_x")+
  scale_x_discrete(drop = TRUE) +
  ylab("Shannon Diversity
       
       ")+
  xlab(" ")

#n asvs
p2<-rich %>%
  ggplot(aes(x=Compartment, y=Observed,  col= Fraction, fill=Fraction))+
  geom_boxplot(alpha=.5) +
  scale_color_manual(values=mycols3) +
  scale_fill_manual(values=mycols3) +
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")+
  facet_wrap(~Fraction, scales = "free_x")+
  scale_x_discrete(drop = TRUE) +
  ylab("Numbers of ASVs")+
  xlab("")

 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")
svg(file="diversity.svg",width = 8, height=7)
#windows(8,7)
require(gridExtra)
grid.arrange(p1, p2, ncol=1)
dev.off()

# summary table
rich %>% group_by(compartment_BCAT) %>% summarise(mean(Shannon))

#####DIVERSITY STATS######
####shannon###
# overall anova
df<-rich %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT, REP)

m1<-lm(Shannon ~ Compartment+ BONCAT + REP+ Compartment*BONCAT,  data = df)
summary(m1)
# no effect of REP. so we dropped REP from the model. 
# rhizo, roots and nod are all different but there is an interaction with BONCAT signal
# we manual subset into comparisons of interest and did a t test. 

# Total DNA
# rhizo vs bulk
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Bulk_Soil", BONCAT=="DNA")

m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)

# viable cell
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots")
m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1) 


# nod vs endo viable
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule"| Compartment=="Roots")
m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1) 


# active
# rhizo vs endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="POS")
m1<-lm(Shannon ~ Compartment,  data = df)
summary(m1)


# roots vs nod
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots"| Compartment=="Nodule", BONCAT=="POS")
m1<-lm(Shannon ~ compartment_BCAT,  data = df)
summary(m1)



#total DNA verse viable
# in rhizo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="BONCAT") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)


#total DNA verse active
# in rhizo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="SYBR") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)


# active verse viable
# in rhizo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="DNA") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)


# in root
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)


# in nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Shannon, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule")
m1<-lm(Shannon ~ BONCAT,  data = df)
summary(m1)



###No. ASVs. ###
# overall anova 
df<-rich %>%
    select(Observed, BONCAT, Compartment, compartment_BCAT, REP)
  m1<-lm(Observed ~ Compartment+ BONCAT + REP+ Compartment*BONCAT,  data = df)
  summary(m1)
# rep was non signifcant so we dropped it from the model.
# rhizo, roots and nod are all different but there is an interaction with BONCAT signal

# Total DNA
# rhizo vs bulk
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Bulk_Soil", BONCAT=="DNA")
  m1<-lm(Observed ~ Compartment,  data = df)
  summary(m1)

# viable cell
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Roots")
  m1<-lm(Observed ~ Compartment,  data = df)
  summary(m1) 

# nod vs endo viable
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT=="SYBR") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule"| Compartment=="Roots")
  m1<-lm(Observed ~ Compartment,  data = df)
  summary(m1) 

# active
# rhizo vs endo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere"| Compartment=="Nodule", BONCAT=="POS")
  m1<-lm(Observed ~ Compartment,  data = df)
  summary(m1)

# roots vs nod
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots"| Compartment=="Nodule", BONCAT=="POS")
  m1<-lm(Observed ~ compartment_BCAT,  data = df)
  summary(m1)


#total DNA verse viable
# in rhizo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="Active_Cell") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
  m1<-lm(Observed ~ BONCAT,  data = df)
  summary(m1)

#total DNA verse active
# in rhizo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="SYBR") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
  m1<-lm(Observed ~ BONCAT,  data = df)
  summary(m1)

# active verse viable
# in rhizo
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", BONCAT!="DNA") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Rhizosphere")
  m1<-lm(Observed ~ BONCAT,  data = df)
  summary(m1)

# in roots
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Roots")
  m1<-lm(Observed ~ BONCAT,  data = df)
  summary(m1)

# in nodule
df<-rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive") %>%
  select(Observed, BONCAT, Compartment, compartment_BCAT)%>%
  filter(Compartment=="Nodule")
  m1<-lm(Observed ~ BONCAT,  data = df)
  summary(m1)

#####import data & rarefy ###########
## beta diversity is done on rarefied reads
## Set the working directory
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/")
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

## Determine minimum available reads per sample ##
min.s<-min(rowSums(asvs.t))
min.s
### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs.t, min.s)
dim(asvs.t)
dim(asvs.r)

###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Viable_Cell', 'POS'='Active_Cell', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Fraction, metadat$Compartment))

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
df<-subset_samples(ps, Compartment=="Roots" | Compartment=="Nodule" )
df<-prune_taxa(taxa_sums(df) > 0, df)

# remove chloroplasts
remove<-subset_taxa(df, Domain=="Unassigned" |  Phyla=="" | Phyla==" p__"  ) 
remove
### These 5 ASVS were not mitochondria or chloroplasts, the remain 48 were
kp<-c("a49a51f3a3e3ea140206b10c5665cc13", "55c30bcbeacfedffa7aeb332600548b2" , "cfeae1df224b7e426ea125ab2bb824fc", "b260024f11a7d77d4f03e5ca2e239860", "f5211207035c7ea5b6f2ecfdad3765e1")
badtaxa<-taxa_names(remove)
badtaxa <- badtaxa[!(badtaxa %in% kp)]
length(badtaxa)

########## remove these guys
alltaxa<-taxa_names(ps)
mytaxa <- alltaxa[!(alltaxa %in% badtaxa)]
ps<-prune_taxa(mytaxa, ps )
ps.r<-prune_taxa(taxa_sums(ps) > 0, ps)
# ps.r is rarefied. ps is not.
# 12652 taxa

# check  reads after QC
n<-rowSums(otu_table(ps.r))
# number of ASVS after QC
n<-otu_table(ps)
n<-rowSums(ifelse(n[]>0,1,0))
n


#####PCOA plots FIG 4######

##all comparents##
#Pcoa on rarefied asvs Data
ps.r<-subset_samples(ps.r, Compartment !="ctl")
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps.r), method = "bray")
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


# subset metadata
metadat2<-filter(metadat, Compartment!="ctl")
as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
unique(levels(as.factor(metadat2$Fraction)))
#color and shapes
mycols_pc <-  c("#fc8449", "#bcd3e8" , "#282c55" )
#get shapes
square <- 22
diamond <- 23
triangle <- 24
circle <- 21
upsidedown_tri <- 25


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")

pdf(file="Pcoa.pdf",width = 8, height=8)
par(mfrow=c(2,2))

#all
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="All Compartments",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
title(adj = 0, main= "A")
points(otus.p, col=c("darkgrey"),
       pch=c(25, circle, diamond, triangle)[as.factor(metadat2$Compartment)],
       lwd=1,cex=2,
       bg=mycols_pc[as.factor(metadat2$Fraction)])
legend("topleft", legend=c( "Active Cell"   ,   "Viable Cell", "TotalDNA"  ),
       fill= mycols3,
       cex=1,
       title = "Fraction",
       bty = "n")
legend("top", legend=c("Nodule", "Root", "Rhizosphere", "Bulk soil"  ),
       pch=c(1, 2, 5, 6),
       cex=1,
       title = "Compartment",     bty = "n")
#dev.off()


##soil##
ps2<-subset_samples(ps.r, Compartment !=  "Nodule" & Compartment != "Roots" & Compartment !="ctl" & Fraction != "Total_DNA")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
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
metadat2<-as.data.frame(sample_data(ps2))
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots" & Compartment!="ctl" & Fraction != "Total_DNA")
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/pcoa/soil_raw.svg",width = 4, height=4 )
#windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Rhizosphere",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
title(adj = 0, main= "B")
points(otus.p[,1:2],
       col=c("darkgray"),
       pch=(diamond),
       lwd=1,cex=2,
       bg= c(mycols3[1:2])[as.factor(metadat2$Fraction)])
ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= mycols3[1:2],
            alpha = 50)

## roots ###
ps2<-subset_samples(ps.r, Compartment == "Roots")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
df<-as.data.frame(otu_table(ps2))
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# Perform PCoA analysis of BC distances #
otus.pcoa <- cmdscale(otus.bray, k=(9-1), eig=TRUE)
# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]
# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]
metadat2<- sample_data(ps2)

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/pcoa/roots.svg",width = 4, height=4 )
#windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Roots",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
title(adj = 0, main= "C")
points(otus.p[,1:2],col=c("darkgrey"),
       pch=triangle,
       lwd=1,cex=2,
       bg=c(fill= c(mycols3[1:2]) )[as.factor(metadat2$Fraction)])
ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= c(mycols3[1:2]),
            alpha = 50)

###nodule###
ps2<-subset_samples(ps.r, Compartment == "Nodule")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
df<-as.data.frame(otu_table(ps2))
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
# subset metadata
metadat2<- as.data.frame(sample_data(ps2))

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/pcoa/nodule.svg",width = 4, height=4 )
#windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Nodule",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
title(adj = 0, main= "D")
points(col=c("darkgrey"),
       otus.p[,1:2],
       pch=c(circle),
       lwd=1,cex=2,
       bg=c(fill= c(mycols3[1:2] ) )[as.factor(metadat2$Fraction)])
ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= c(mycols3[1:2]),
            alpha = 50)

dev.off()

#####PCOA STATS----------------#########

###BETA DISPERSION#
#full data set between compartments#
asvs.clean<-otu_table(ps.r)
# subset metadata
metadat2<-filter(metadat, Compartment!="ctl")
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(otu_table(ps.r), method = "bray")

dispersion <- betadisper(asvs.bray, group=metadat2$Compartment)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     3 0.037320 0.0124400 22.441    999  0.001 ***

#full data between fractions
dispersion <- betadisper(asvs.bray, group=metadat2$Fraction)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#Response: Distances
#          Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.01847 0.009235 0.1262    999  0.878
#Residuals 37 2.70752 0.073176   

# just soil
ps2<-subset_samples(ps.r, Compartment !=  "Nodule" & Compartment != "Roots" & Compartment !="ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots" & Compartment!="ctl") 
#test for dispersion between compartments
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$compartment_BCAT))
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     3 0.037320 0.0124400 22.441    999  0.001 ***


# just rhizosphere
ps2<-subset_samples(ps.r, Compartment ==  "Rhizosphere" & Compartment !="ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment ==  "Rhizosphere" & Compartment!="ctl") 
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
anova(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     2 0.033450 0.0167251 29.276    999  0.001 ***
#  Residuals 11 0.006284 0.0005713   





# just rhizosphere total cells and active
ps2<-subset_samples(ps.r, Compartment ==  "Rhizosphere" & Compartment !="ctl" & Fraction!="Total_DNA")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)

# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment ==  "Rhizosphere" & Compartment!="ctl" & Fraction!="Total_DNA")
#test for dispersion between groups
# compartments are not equaly dispersed
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999

#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
#Groups     1 0.013346 0.0133460 19.805    999  0.001 ***
#  Residuals  7 0.004717 0.0006739              


# just Total DNA
ps.r
metadat
ps2<-subset_samples(ps.r, Fraction=="Total_DNA" & Plant =="Clover")
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
#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#Response: Distances
#Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.00092857 0.00092857 2.3009    999  0.166
#Residuals  7 0.00282501 0.00040357 
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse



#roots + nodules #
ps2<-subset_samples(ps.r, Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
otus.bray<-vegdist(otu_table(ps2), method = "bray")
metadat2<-metadat%>% filter(Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse



#roots + nodules #
ps2<-subset_samples(ps.r, Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
otus.bray<-vegdist(otu_table(ps2), method = "bray")
metadat2<-metadat%>% filter(Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Compartment))
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse



#roots #
ps2<-subset_samples(ps.r, Compartment ==  "Roots" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
otus.bray<-vegdist(otu_table(ps2), method = "bray")
metadat2<-metadat%>% filter(Compartment == "Roots")
metadat2
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse


#nodules #
ps2<-subset_samples(ps.r, Compartment ==  "Nodule" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
otus.bray<-vegdist(otu_table(ps2), method = "bray")
metadat2<-metadat%>% filter(Compartment == "Nodule")
#test for dispersion between groups
dispersion <- betadisper(otus.bray, group=as.factor(metadat2$Fraction))
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE) ##sd ellipse

###PERMANOVA##
#full model
asvs.clean<-otu_table(ps.r)
dim(asvs.clean)
head(asvs.clean)
head(metadat)
metadat2<-filter(metadat, Compartment!="ctl")
# Calculate Bray-Curtis distance between samples
asvs.bray<-vegdist(otu_table(ps.r), method = "bray")
#
asvs.perm<- adonis2(asvs.clean ~ Compartment+ Fraction+Compartment*Fraction +REP, data = metadat2, permutations = 999, method="bray")
asvs.perm
#Df SumOfSqs      R2       F Pr(>F)    
#Compartment           3   7.9377 0.68157 33.0965  0.001 ***
#  Fraction              2   0.8828 0.07580  5.5211  0.001 ***
#  Compartment:Fraction  2   0.2675 0.02297  1.6728  0.122    
# rep was not sifnificant so we dropped it from the model.

#rhizo +bulk total dna
ps2<-subset_samples(ps.r, Compartment ==  "Rhizosphere" | Compartment =="Bulk_Soil" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment =="Rhizosphere" | Compartment =="Bulk_Soil" )
otu.perm<- adonis2((otu_table(ps2))~ Compartment, data = metadat2, permutations = 999, method="bray")
otu.perm


# root + nodules 
ps2<-subset_samples(ps.r, Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment !=  "Rhizosphere" & Compartment !="Bulk_Soil" & Compartment != "ctl")
#full model
otu.perm<- adonis2((otu_table(ps2))~ Compartment+ Fraction + Compartment*Fraction, data = metadat2, permutations = 999, method="bray")
otu.perm

# roots + nodules, difference between compartments
otu.perm<- adonis2((otu_table(ps2))~ Compartment, data = metadat2, permutations = 999, method="bray")
otu.perm

# roots + nodules, difference between fractions
otu.perm<- adonis2((otu_table(ps2))~ Fraction, data = metadat2, permutations = 999, method="bray")
otu.perm


# roots difference between fraction
ps2<-subset_samples(ps.r, Compartment =="Roots" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment =="Roots")
otu.perm<- adonis2((otu_table(ps2))~ Fraction , data = metadat2, permutations = 999, method="bray")
otu.perm
#adonis2(formula = (otu_table(ps2)) ~ Fraction, data = metadat2, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)  
#Fraction  1  0.14573 0.53911 8.1879  0.011 *


# nodule difference between fractions
# roots difference between fraction
ps2<-subset_samples(ps.r, Compartment =="Nodule" )
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
# Calculate Bray-Curtis distance between samples
otus.bray<-vegdist(otu_table(ps2), method = "bray")
# subset metadata
metadat2<-metadat%>% filter(Compartment =="Nodule")
otu.perm<- adonis2((otu_table(ps2))~ Fraction , data = metadat2, permutations = 999, method="bray")
otu.perm
#adonis2(formula = (otu_table(ps2)) ~ Fraction, data = metadat2, permutations = 999, method = "bray")
#Df SumOfSqs      R2      F Pr(>F)  
#Fraction  1 0.065957 0.51041 7.2978  0.031 *


#####SCATTERPLOT FIG 5##############
# this is done on rarefied data, to be able to compare abundances between samples
sample_data(ps.r)
df<-subset_samples(ps.r, Compartment!="ctl" & Compartment=="Rhizosphere" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df<-as.data.frame(t(otu_table(df)))
dim(df)
# 5215 taxa

colnames(df)
n<-c( "BONCAT_1", "Total_1", "BONCAT_2" ,  "Total_2",  "BONCAT_3" , 
      "Total_3" , "BONCAT_4" ,  "Total_4",  "Total_5" )
colnames(df)<-n
#make a column for the value in active and a column for the value in viable
#make a columns that states the rep. 
df$otu <- row.names(df) 
total<-df %>% dplyr::select(contains("Total"))
total$otu <- row.names(total) 
total<-total %>% pivot_longer(cols = 1:4, values_to = "total", names_to = "rep_total") %>% dplyr::select(-Total_5)
#active
active<-df %>% dplyr::select(contains("BONCAT"))
active$otu <- row.names(active) 
active<-active %>% pivot_longer(cols = 1:4, values_to = "active", names_to = "rep")
active[,2:3]
df<-cbind(total, active[,2:3])
df<-df%>% filter(total>0 & active>0 )
df<-df%>%
  mutate( rep= str_split_i(df$rep, "_", 2))

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig5_scatter")
mycols<-c("#00436d", "#416d9a", "#729ac9", "#a3c9fa")

svg(file="scatterplot.svg",width = 7, height=4)
windows(8,4)

filter(df) %>%
  ggplot(aes(x=log10(total) , y=log10(active), col=rep)) + 
  geom_jitter(stat="identity", position="identity") +
  #geom_abline(slope=1)+
  geom_smooth(method=lm, se=FALSE)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = mycols)+
  labs(x="log10 Abundance in Total Viable Cell population",y="log10 Abundance in Active population")
dev.off()
#lm
m1<-lm(log10(total)~log10(active), data=df)
summary(m1)
#R squared = .32


#####BARPLOT TOP ASVS FIG 5######
# this is on rarefied ASVS
# filter for total and active

df<-subset_samples(ps.r, Compartment!="ctl"& Compartment!="Bulksoil" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df<-as.data.frame(t(otu_table(df)))
df<-cbind(taxon,df)

#get means
df$nodule.total.mean <-   rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
df$nodule.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root.total.mean <-   rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
df$root.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
colnames(df)
#remove extra columns
asv <- dplyr::select(df, -contains("POS") )  %>% dplyr::select(., -contains("SYBR"))
head(asv)

# select top 50 ASVs viable rhizosphere
asv<-asv[order(asv$rhizo.total.mean, decreasing = TRUE),]
top<-asv[1:50,]
top$asv<-row.names(top)
top

# what percent of the whole community are the top asvs
all<-colSums(asv[,12:13])
t<-colSums(top[,12:13])
t/all
# top 50 otus is 27% of the total population;

#edit labels
top$Phyla<-sub("p__", "", top$Phyla)
top<-top[order(top$rhizo.total.mean, decreasing = TRUE),]
top$otu1<-c(1:50)
top$otu1<-paste0("ASVS ",top$otu1)
top

#colors
# acido - dk blue
# actino - light blue
# bacteriodota - pink
# chlorofexi light red
# proteobacteria, light green
# Methylomirabilota
# verrucomicrobiota - pale gold


mycols8<- c( "#1F78B4","#A6CEE3", "#75026d",  "#ed6361",  "#6A3D9A", "#B2DF8A", "#FF7F00","#FDBF6F")

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/")

svg(file="Fig5_scatter/barplot_toptaxa.svg",width = 8, height=5)
windows(6,5)
ggplot(top)+
  geom_bar(aes(x=reorder(otu1, -rhizo.total.mean) , y=rhizo.bcat.mean, fill= Phyla), stat="identity", position="identity") +
  geom_line(aes(x=as.numeric(reorder(otu1, -rhizo.total.mean)) , y=rhizo.total.mean), linewidth=1) +
  xlab("top 50 Viable Cell ASVS")+
  ylab("Average rarefied reads per sample")+
  scale_fill_manual(values=mycols8)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_blank(), legend.position = c(.8, .7))+
  guides(fill=guide_legend(title="Abundance Active Cells"))
dev.off()


## with labs
svg(file="Fig5_scatter/barplot_toptaxa_lab.svg",width = 8, height=5)
ggplot(top)+
  geom_bar(aes(x=reorder(otu1, -rhizo.total.mean) , y=rhizo.bcat.mean, fill= Phyla), stat="identity", position="identity") +
  geom_line(aes(x=as.numeric(reorder(otu1, -rhizo.total.mean)) , y=rhizo.total.mean), linewidth=1) +
  xlab("top 50 Viable Cell ASVS")+
  ylab("Average rarefied reads per sample")+
  scale_fill_manual(values=mycols8)+
  theme_classic(base_size = 14)+
  theme(axis.text.x=element_text(angle=90), legend.position = c(.8, .7))+
  guides(fill=guide_legend(title="Abundance Active Cells"))
dev.off()
    
#####BARPLOT PHYLA LEVEL FIG 6 #######
#use non rarefied taxa because we are analyzing the proportion of read we will normalize by Number of reads.
# filter for total and active rhizosphere
    ps
    sample_data(ps)
    df<-subset_samples(ps, Compartment!="ctl"& Compartment=="Rhizosphere" & Fraction!="Total_DNA")
    df<-prune_taxa(taxa_sums(df) > 0, df)
    taxon<-tax_table(df)
    df1<- as.data.frame(otu_table(df)) # this is for later :)
    df<-as.data.frame(otu_table(df))
    dim(df)
    
# normalize by number of reads
    df<-df/rowSums(df)
    df<-as.data.frame(t(df))
    head(df)
    # 4854 taxa
    df<-cbind(taxon,df)
 
    #aggregate
    df<-aggregate(cbind(C10R.POS_S30, C10R.SYBR_S20, C1R.POS_S27,   C1R.SYBR_S16, C2R.POS_S28 , 
           C2R.SYBR_S17 , C5R.POS_S29 ,  C5R.SYBR_S18 , C7R.SYBR_S19  ) ~ Phyla, data = df, FUN = sum, na.rm = TRUE)
    colnames(df)<- c("Phyla", "BONCAT_1" , "Total_1", "BONCAT_2",   "Total_2" , "BONCAT_3"  , "Total_3" ,
                    "BONCAT_4"  , "Total_4",  "Total_5") 
  
    head(df)
    dim(df)
    
# grab top phyla
    row.names(df) <- df$Phyla
    df$sum <- rowSums(df[,2:9])
    df<-df[order(-df$sum),]    
    top<-df[c(1:5,7),]
  # rm sum
  top
  top<-top%>%
  pivot_longer(2:10, values_to = "abundance", names_to= "rep_fraction" )
  top<-top%>%mutate(fraction = str_split_i(top$rep_fraction, "_", 1)) %>%
  mutate(Phyla= str_split_i(top$Phyla, "_", 3))
  head(top)
  
#label
top<-top %>% mutate(fraction=recode(fraction, 'Total'='Viable Cells'))
top$fraction<-factor(top$fraction, levels = c('Viable Cells', 'BONCAT'))
#colors
# acido - dk blue
# actino - light blue
# bacteriodota - pink
# chlorofexi light red
# proteobacteria, light green
# Methylomirabilota
# verrucomicrobiota - pale gold


mycols8<- c( "#1F78B4","#A6CEE3", "#75026d",  "#ed6361",  "#6A3D9A", "#B2DF8A", "#FF7F00","#FDBF6F")

phycols7<-c("#1F78B4","#A6CEE3", "#75026d",  "#FB9A99", "#33A02C","#FF7F00",  "#FB9A99")

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")
svg(file="barplot_phyla.svg",width = 8, height=4)
windows(7,3)
ggplot(top) +
  geom_boxplot(aes(x= reorder(fraction, -sum), y= abundance, fill=Phyla), outlier.shape = NA, alpha=.5, show.legend = FALSE)+
  geom_jitter(aes(x= reorder(fraction, -sum), y= abundance, fill=Phyla), position = position_jitter(width = .2), show.legend = FALSE)+
  facet_wrap(~reorder(Phyla, -sum), nrow=1)+
  theme_minimal(base_size = 12)+
  scale_fill_manual(values=phycols7)+
  labs(y="proportion of reads in sample")+
  theme(axis.text.x=element_text(angle=45, hjust=0.9), axis.title.x = element_blank())
dev.off()
#####BARPLOT STATS ########
# mixed model with binomial regression
# we need a data frame of # success to # failure


df<-subset_samples(ps, Compartment!="ctl"& Compartment=="Rhizosphere" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df1<- as.data.frame(otu_table(df)) # this is for later :)
# make column that is the number of reads
reads<-rowSums(df1)
df1<-as.data.frame(t(df1))
# add taxa info
df1<-cbind(taxon,df1)
#aggregate
df1<-aggregate(cbind(C10R.POS_S30, C10R.SYBR_S20, C1R.POS_S27,   C1R.SYBR_S16, C2R.POS_S28 , 
                    C2R.SYBR_S17 , C5R.POS_S29 ,  C5R.SYBR_S18 , C7R.SYBR_S19  ) ~ Phyla, data = df1, FUN = sum, na.rm = TRUE)

# grab top phyla
row.names(df1) <- df1$tPhyla
df1$sum <- rowSums(df1[,2:9])
df1<-df1[order(-df1$sum),]    
top<-df1[c(1:5,6,7),]
top
# pivot
top<-top%>%
  pivot_longer(2:10, values_to = "abundance", names_to= "rep_fraction" )
top<-top%>%mutate(fraction = str_split_i(top$rep_fraction, "_", 1)) %>%
  mutate(Phyla= str_split_i(top$Phyla, "_", 3))

#put No. reads in
top$reads =rep(reads, 7 )
top<-select(top , -sum, -fraction, -rep_fraction)
top$failures = top$reads-top$abundance
trt<-c(rep(c("Active", "Total"), 4), "Total")
top$trt = rep(trt, 7)
top$rep<-rep(c(1,1, 2, 2, 3, 3,4, 4,5),7)
top


#proteobacteria model
df<-filter(top, Phyla=="Proteobacteria") %>% select(abundance,failures, trt, rep )

m1<-glm(data= df, cbind(abundance,failures)~trt+rep, family = binomial)
anova(m1, test= "LRT")  
#
#acido model
df<-filter(top, Phyla=="Acidobacteriota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df, cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  


#Verrucomicrobiota model
df<-filter(top, Phyla=="Verrucomicrobiota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  
#not sig

# Actinobacteriota model
df<-filter(top, Phyla=="Actinobacteriota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  


#Bacteroidota model
df<-filter(top, Phyla=="Bacteroidota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  
#not sig 

#Planctomycetota model
df<-filter(top, Phyla=="Planctomycetota") %>% select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
m1
anova(m1, test= "LRT")  

#Chloroflexota model
df<-filter(top, Phyla=="Chloroflexota") %>%select(abundance,failures, trt, rep )
m1<-glm(data= df,cbind(abundance,failures)~trt+rep, family = binomial)
anova(m1, test= "LRT")  



##### ANCOM ###############
# ancom is run on unrarefied data
# subset to rhizosphere active and viable
ps1<-subset_samples(ps, Compartment=="Rhizosphere" & BONCAT!="DNA")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# add 1 to everything... b/c ancom can't deal with structural zeros.
df<-as.data.frame(otu_table(ps1))
df<-df+1
# remove other column
df<-select(df,-Others)

# add asv column
taxon<-as.data.frame(tax_table(ps1))
taxon$asv<-row.names(taxon)
#remove others
taxon<-filter(taxon, asv!="Others")

Workshop_OTU <- otu_table(as.matrix(df), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(ps1)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps1 <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps1

### skip to data import if you already ran ancom
sample_data(ps1)

ps2 = mia::makeTreeSummarizedExperimentFromPhyloseq(ps1)

ps2

out1 = ancombc(data = ps2, assay_name = "counts", 
               tax_level = "asv", phyloseq = NULL, 
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
               group = "Fraction", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
sample_data(ps1)

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/tables")
#log fold change
tab_lfc = res$lfc
head(tab_lfc)
dim(tab_lfc)
col_name = c("asv", "LFC_Intercept", "LFC_FractionViable_Cell")
colnames(tab_lfc) = col_name
head(tab_lfc)
write_delim(as.data.frame(tab_lfc), file = "ancom_Log_fold_change.txt", delim = " ")

# standard error
tab_se = res$se
col_name = c("asv", "se_Intercept", "se_viable")
colnames(tab_se) = col_name
head(tab_se)
tab_se<-as.data.frame(tab_se)
write.table(tab_se, file = "ancom_SE.txt")
tab_se<-read.table("ancom_SE.txt", header = TRUE)
head(tab_se)

#test statistcs W maybe it's willcoxin?
tab_w = res$W
col_name = c("asv", "W_Intercept", "W_viable")
colnames(tab_w) = col_name
head(tab_w)
write.table(as.data.frame(tab_se), file = "ancom_SE.txt")

# P-values from the Primary Result
tab_p = res$p_val
col_name = c("asv", "p_Intercept", "p_viable")
colnames(tab_p) = col_name
head(tab_p)
write.table(as.data.frame(tab_p), file = "ancom_pval.txt")


#Adjusted p-values from the Primary Result"
tab_q = res$q
head(tab_q)
col_name = c("asv", "adj_p_Intercept", "adj_p_viable")
colnames(tab_q) = col_name
head(tab_q)
write.table(as.data.frame(tab_se), file = "ancom_adjpval.txt")

# yes or no is a taxa differentially abundant
tab_diff = res$diff_abn
col_name = c("asv", "DA_Intercept", "DA_Fraction_Total_cells_Active")
colnames(tab_diff) = col_name
head(tab_diff)
write_delim(as.data.frame(tab_diff), file = "ancom_DA.txt", delim = " ")

# through all togetha nd remove anything with DNA
tab <-tab_lfc %>%
    left_join(., tab_se) %>%
    left_join(., tab_w ) %>%
    left_join(., tab_p) %>%
    left_join(., tab_q) %>%
    left_join(., tab_diff)
write.table(as.data.frame(tab), file = "ancom_table.txt")

######import df from ANCOM######
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/tables")
tab<-read.table("ancom_table.txt", header = TRUE)
head(tab)

##add add abundance and taxon info
df<-as.data.frame(t(otu_table(ps1)))
df$asv<-row.names(df)
df<-left_join(df, tab)
taxon<-as.data.frame(tax_table(ps1))
df<-left_join(taxon, df)
df

### summarise the abundance in active and total
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
t<-df %>% select(contains("R.SYBR"))
sd_total<- apply(t, 1, sd, na.rm=TRUE)
sd_total
df$sd_total <- sd_total
  
df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
t<-df %>% select(contains("R.POS"))
sd_active<- apply(t, 1, sd, na.rm=TRUE)
sd_active
df$sd_active <- sd_active
head(df)

#remove some columns to make things simpler
df<-df %>% select(-contains("R."))
df<-df %>% select(-contains("Intercept"))
df<-df %>% select(-contains("DNA"))
colnames(df)

#quick histogream
hist(df$rhizo.bcat.mean, breaks = 100)
hist(df$rhizo.total.mean, breaks=100)

# make a column the lowest taxa rank assigned 
df$label<-df$Species
df$label[which(df$label==" s__"  )] <- df$Genus[which(df$label==" s__" )]
df$label[which(df$label==""  )] <- df$Genus[which(df$label=="" )]
df$label[which(df$label==" g__"  )] <- df$Family[which(df$label==" g__" )]
df$label[which(df$label==""  )] <- df$Family[which(df$label=="" )]
df$label[which(df$label==" f__"  )] <- df$Order[which(df$label==" f__" )]
df$label[which(df$label== " g__SCN-69-37" )] <- df$Order[which(df$label==" g__SCN-69-37" )]
df$label[which(df$label== " f__UBA2999" )] <- df$Order[which(df$label== " f__UBA2999"  )]
df$label[which(df$label== " g__WHSN01" )] <- df$Order[which(df$label== " g__WHSN01"  )]
df$label[which(df$label== " s__UBA11740 sp003168335"    )] <- df$Order[which(df$label== " s__UBA11740 sp003168335"    )]
df$label[which(df$label== " g__PSRF01"      )] <- df$Order[which(df$label== " g__PSRF01"    )]
df$label[which(df$label==    " s__VFJQ01 sp009885995"    )] <- df$Family[which(df$label==  " s__VFJQ01 sp009885995"  )]
df$label[which(df$label==   " s__OLB17 sp001567505"       )] <- df$Family[which(df$label==  " s__OLB17 sp001567505"     )]

# right now the active community is the reference. It's a little confusing. 
# multiply by -1 to make it so the viable community is the reference! 
# so a negative number would be depleted in active
df$LFC_FractionViable_Cell<-df$LFC_FractionViable_Cell*-1

#rm "p__" in phyla
df$Phyla<-sub(" p__", "", df$Phyla)
# number ASvs
df$asv_no <- paste0("ASV", row.names(df))

# quick overall volcano plot to check distribution
ggplot(df, aes(x=LFC_FractionViable_Cell , y=p_viable)) + 
  geom_jitter()+ 
  theme_bw()
# lfc verse the pvalue
ggplot(df, aes(x=LFC_FractionViable_Cell , y=-log(p_viable), col=DA_Fraction_Total_cells_Active)) + 
  geom_jitter()+ 
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  theme_bw( )

#LOOK AT DA taxa
DA<-df %>% filter(DA_Fraction_Total_cells_Active==TRUE)
dim(DA)
head(DA)
#save DA taxa file
write.csv(DA, file="ancom_DA_taxa.csv")

# sort for negative slope - more dormant taxa
dormant<-DA[DA$LFC_FractionViable_Cell <0,]
dormant<-dormant[order(dormant$rhizo.total.mean, decreasing = TRUE),]
dormant
write.csv(dormant, file ="ancom_dormanttaxa.csv")

# sort for more active taxa - postive slope
active<-DA[DA$LFC_FractionViable_Cell > 0,]
active<-active[order(active$rhizo.bcat.mean, decreasing = TRUE),]
head(active)
write.csv(active, file ="ancom_activetaxa.csv")

######DA TAXA figure ############ 
# filter for at least 50 reads
DA<-filter(DA, DA$rhizo.bcat.mean>50 |  DA$rhizo.total.mean>50  )
dim(DA)

#remove numbers as the end of labels to make them easier to read in the figure
DA$label<-gsub("_48326", "", DA$label)
DA$label<-gsub("_48670", "", DA$label)
DA$label<-gsub("_A_50105", "", DA$label)

#add asvs No. label
DA$label<-paste(DA$asv_no, DA$label)
DA$label

DA$type <- "Differentially Abundant ASVS"

#add * for DA taxa
DA$sig<-ifelse( DA$DA_Fraction_Total_cells_Active=="TRUE", "*", "")
DA$label<-paste(DA$label, DA$sig)
DA$label

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")

### cols for this plot
# actinobacteria - 
# acidobacteria - dk blue 
# actino - light blue
# proteobaceria - light green
# gemma - 


svg(file="barplot_top20_DA.svg",width = 8, height=7.5)
windows(8,7.5)
ggplot(DA)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla),
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=LFC_FractionViable_Cell-se_viable,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values=c("#A6CEE3", "#75026d","#B2DF8A"))+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab("Differentially Abundant Asvs")+
  ylab("Log fold change viable to active")
  #theme(legend.position = "none")
  
dev.off()  


###LFC top 10 taxa ##
## filter for top 50 most abundant taxa or DA taxa
df<-df %>%
  arrange( -rhizo.total.mean)
top<-df[1:10,]

#remove the " f__ " part
top$label <-substr(top$label, 5, nchar(top$label)-1)

#remove numbers as the end of labels to make them easier to read in the figure
top$label<-gsub("_48326", "", top$label)
top$label<-gsub("_48670", "", top$label)
top$label<-gsub("_A_50105", "", top$label)
top$label<-gsub("_E_64746", "", top$label)
top$label<-gsub("_58024", "", top$label)
top$label<-gsub("_A_58049", "", top$label)
top$label

#add asvs label
top$label<-paste(top$asv_no, top$label)

#add * for DA taxa
top$sig<-ifelse( top$DA_Fraction_Total_cells_Active=="TRUE", "*", "")
top$label<-paste(top$label, top$sig)
top$label

#group label
top$type <- "Most Abundant Asvs"
### cols for this plots
# acidobacteria - dk blue 
# proteobaceria - light green
head(top)
#plot

svg(file="barplot_top10.svg",width = 8, height=3.5)
windows(8,8)
ggplot(top)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla), 
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=LFC_FractionViable_Cell-se_viable,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values=c("#1F78B4", "#B2DF8A"))+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab("Most abudant Asvs")+
  ylab("Log fold change viable to active")+
  # facet_grid( scales = "free", space = "free",  rows=vars(Phyla)) 
  theme(legend.position = "none")
dev.off()  

##DA fig combined###
#combine
combine<-full_join(top, DA)
head(combine)

#colors
"#06568c" #  acidobacteria - dk blue 
"#B2DF8A" # proteobaceria - light green
"#52b8d1" # actinobacteria - light blue
"#d40d63" # Gemmatimonadota raspberry

mycols<-c("#06568c",   "#52b8d1",   "#d40d63", "#B2DF8A" )
#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")
svg(file="barplot_big.svg",width = 8, height=7)
windows(8,8)
ggplot(combine)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla), 
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values= mycols)+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab(" Asvs")+
  ylab("Log fold change viable to active")
  #facet_grid( space = "free",  rows=vars(Type))
  #theme(legend.position = "none")
dev.off()  

#plot
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig6_DAtaxa")
svg(file="barplot_big_label.svg",width = 8, height=7)
windows(8,8)
ggplot(combine)+
  geom_bar(aes(x=reorder(label, +rhizo.total.mean) , y=LFC_FractionViable_Cell, fill= Phyla), 
           stat="identity", position="dodge")+
  geom_errorbar(aes(x=label, ymin=-se_viable+LFC_FractionViable_Cell,
                    ymax=LFC_FractionViable_Cell+se_viable))+ 
  scale_fill_manual(values= mycols)+
  theme_minimal(base_size = 14) +
  coord_flip()+
  ylim(c(-5,3))+
  xlab(" Asvs")+
  ylab("Log fold change viable to active")+
  #facet_grid( space = "free",  rows=vars(Type))+
  theme(legend.position = "none")
dev.off()  

#####VENN DIAGRAM Viable#####
#on rarefied data
sample_data(ps.r)
ps1<-subset_samples(ps.r, Fraction=="Viable_Cell")
# at least in 3 samples min reads is 50
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 542 taxa

# grab taxonomy + asv column
taxon<-as.data.frame(tax_table(ps1))
taxon$asv<-row.names(taxon)
taxon<-filter(taxon, asv!="Others")

#grab data and rename rows by compartment
df<-as.data.frame(otu_table(ps1))

head(df)
n<-row.names(df)
n[grepl("N" , n)]="nodule"
n[grepl("E" , n)]="roots"
n[grepl("R" , n)]="rhizo"
n

# sum by compartment
df<-rowsum(df, n)
#transform
df<-as.data.frame(t(df))
#remove others row
df<-df[-which(row.names(df)=="Others"),]
head(df)

# present = 1
df[df>1] <- 1

# make 3 groups
nodule<-rownames(df[df$nodule==1,])
roots<-rownames(df[df$roots==1,])
rhizo<-rownames(df[df$rhizo==1,])

x <- list(
  nodule = nodule, 
  roots = roots, 
  rhizo = rhizo
)


#### venn diagram #
library(ggvenn)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig7_heatmap")
svg(file="OTU_level_total_venn.svg",width = 4, height=4 )
windows(4,4)
mycols= c( "#a3c9fa", "#4e97ed", "#045bc2")
ggvenn(
  x, 
  fill_color = mycols,
  stroke_size = 2, set_name_size = 4, text_size = 3, digits = 1, fill_alpha=.6
  #auto_scale = TRUE
) +
  ggtitle("Viable Cells")
dev.off()



#####VENN DIAGRAM Active######
### need 50 reads + in 3 samples
sample_data(ps)
ps1<-subset_samples(ps, Fraction=="Active_Cell")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
# 324 taxa

# grab taxonomy + asv column
taxon<-as.data.frame(tax_table(ps1))
taxon$asv<-row.names(taxon)
taxon<-filter(taxon, asv!="Others")

#grab data and rename rows by compartment
df<-as.data.frame((otu_table(ps1)))
head(df)
n<-row.names(df)
n[grepl("N" , n)]="nodule"
n[grepl("E" , n)]="roots"
n[grepl("R" , n)]="rhizo"
n

# summ by compartment
df<-rowsum(df, n)
df<-as.data.frame(t(df))
head(df)
#remove others
df<-df[-which(row.names(df)=="Others"),]

#make present =1 and split into groups
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

#######venn diagram #
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig7_heatmap")
library(ggvenn)

svg(file="active_venn.svg",width = 4, height=4 )

mycols = c( "#e89c6f", "#d1663f", "#992600")

windows(4,4)
ggvenn(
  x, 
  fill_color = mycols,
  stroke_size = 2, set_name_size = 4, text_size = 3, digits = 1, fill_alpha=.7
  #auto_scale = TRUE
 ) +
  ggtitle("Active Otus")
dev.off()



#####HEATMAP#########
##big heatmap, no phylogeny ##
# df of values in total and active
# at least in 3 samples min reads is 50
sample_data(ps.r)
df<-subset_samples(ps.r, Fraction=="Active_Cell" | Fraction=="Viable_Cell")
df<-ps_prune(df, min.samples = 3, min.reads = 50)
df # 1018 taxa
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
head(df)

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
df$Family[grepl(" f__Bacillaceae_H", df$Family)] <- " Bacillaceae_H"
df$Family[grepl(" f__Streptomycetaceae", df$Family)] <-  " f__Streptomycetaceae"

#aggregate
df<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,  C1R.SYBR_S16,  C2E.SYBR_S22,
                    C2N.SYBR_S15,  C2R.SYBR_S17 ,  C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25,  C7R.SYBR_S19,
                    C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32, C2N.POS_S62,  C5E.POS_S33,
                    C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 156 families

### summarize by compartment 
df$Viable_rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB"))) %>%   glimpse()
df$Viable_root <-   rowMeans(df %>% dplyr::select(contains("E.SYB"))) %>%   glimpse()
df$Viable_nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB"))) %>%   glimpse()

#### summarize by compartment 
df$Active_rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$Active_root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df$Active_nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()

df<-df%>% select(c(Active_rhizo, Active_root, Active_nodule, Viable_rhizo,Viable_root, Viable_nodule )) %>% mutate(otu= row.names(df))
head(df)

#rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>%   glimpse()
df<-mutate(df, Family=otu) %>% select(., -otu)
head(df)
#rm families that are unknown

df<-df[ !grepl(" f__UBA", df$Family) , ]
df<-df[ !grepl(" f__SG8", df$Family) , ]
df<-df[ !grepl(" f__SCT", df$Family) , ]
df<-df[ !grepl(" f__AC-14", df$Family) , ]
df<-df[ !grepl(" f__SCN", df$Family) , ]
df<-df[ !grepl(" f__RSA", df$Family) , ]
df<-df[ !grepl(" f__WHT", df$Family) , ]
df<-df[ !grepl(" f__RBG", df$Family) , ]
df<-df[ !grepl(" f__CSP", df$Family) , ]
df<-df[ !grepl(" f__DSM", df$Family) , ]
df<-df[ !grepl(" f__TK", df$Family) , ]
df<-df[ !grepl(" f__QHB", df$Family) , ]
df<-df[ !grepl(" f__2013", df$Family) , ]
df<-df[ !grepl(" f__B-17", df$Family) , ]
df<-df[ !grepl(" f__JA", df$Family) , ]
df<-df[ !grepl(" f__J0", df$Family) , ]
df<-df[ !grepl(" f__Gp", df$Family) , ]
df<-df[ !grepl(" f__GWC", df$Family) , ]
df<-df[ !grepl(" f__Fen", df$Family) , ]
df<-df[ !grepl(" f__FW", df$Family) , ]
df<-df[ !grepl(" f__HR", df$Family) , ]



m <- df
# remove those "F__
m$Family<-sub(" f__", "", m$Family)

row.names(m) <- m$Family

#rm family column
m<-select(m, -Family)

#make matrix

m<-log10(m)
m[m== "-Inf"] <- 0
m
#m<-m[,c(3,2,1)]
m<-as.matrix(m)
m

# clustering agrorithm
distance = dist(m, method = "euclidean")
distance
cluster = hclust(distance, method = "ward.D2")
my_list<-cluster$labels
#### use this order for heatmaps with out phylogeny
m<-m[match(row.names(m),as.character(my_list)),]

#p<-ggtree(cluster) + 
#  geom_tiplab(size=3, align=TRUE, linesize=0, offset = -.2) + 
#  theme_tree2()+
#  xlim_tree(1) 

#gheatmap(p, m, 
#         colnames=FALSE,
#         legend_title="active taxa", offset = .5) 
  #scale_x_ggtree() + 
  #scale_fill_gradient(low= "#fffaa2", high =  "#bb0000", aesthetics = "fill", na.value = "white",
  #                    name="Abundance in Active")+

   #ggtitle("Families in Active Community")
#dev.off()


# creates a own color palette
my_palette <- colorRampPalette(c("white", "#95cefc", "#04063b"))(n = 99)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.01,length=2),  # for red
               seq(0.1,0.8,length=48),           # for yellow
               seq(0.81,5.1,length=50))             # for green

windows(7,7)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")

svg(filename = "big_blue_heatmap.svg", width = 9, height = 12)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = 1,
         lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
         lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
         lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
        #key.title = "log 10 Abundance")
        #key.xlab = "")
         key = FALSE) 
#title("Active", line= -4)
dev.off()  
#key
svg(filename = "key.svg", width = 10, height = 10)

windows(10,10)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          #labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = .00001,
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
          lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
          key.title = "log 10 Abundance")
dev.off()  
#key = FALSE

#### hehe make a red one too
# creates a own color palette
my_palette <- colorRampPalette(c("white", "#fffaa2", "#bb0000"))(n = 99)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.01,length=2),  # for red
               seq(0.1,0.8,length=48),           # for yellow
               seq(0.81,5.1,length=50))             # for green

windows(7,7)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")

svg(filename = "big_red_heatmap.svg", width = 9, height = 12)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = 1,
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
          lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
          #key.title = "log 10 Abundance")
          #key.xlab = "")
          key = FALSE) 
#title("Active", line= -4)
dev.off()  

#key
svg(filename = "red_key.svg", width = 8, height = 8)

windows(10,10)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "none",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          #labRow = cluster$labels,
          margins = c(0,0),
          cexRow= 1,
          cexCol = .00001,
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ) , 
          lhei = c(0.43, 2.6, 0.6) , # Alter dimensions of display array cell heighs
          lwid = c(0.6, 4, 0.6) , # Alter dimensions of display array cell widths
          key.title = "log 10 Abundance")
dev.off()  
          #key = FALSE) 

#####HEATMAP with phylogeny:SUPPLEMENT######
##select active taxa 
# at least in 3 samples min reads is 50
df<-subset_samples(ps.r, Fraction=="Active_Cell")
df<-ps_prune(df, min.samples = 3, min.reads = 50)
df # 510 taxa
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
#remove "other" ASVS where rare taxa went
df<-filter(df, !grepl("Others", row.names(df)))
# remove "other" row where rare taxa went
taxon<-taxon[c(which(row.names(taxon)!="Others")),]

####### agregate to the family level
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)
head(df)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"
df$Family[grepl(" f__Blastocatellaceae", df$Family)] <- " f__Blastocatellaceae"

#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#aggregate
df<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 86 families

#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)

#rm unknown families
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")
df<-filter(df, !grepl("UBA", df$otu))
df<-filter(df, !grepl("SG", df$otu))
df<-filter(df, !grepl("SCN", df$otu))

# rename column
df<-mutate(df, Family=otu) %>% select(., -otu)

#removed 'f__'
df$Family<-sub(" f__", "", df$Family)
row.names(df) <- df$Family

#importtree
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")

#match tree to my taxa
length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
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
dim(df)
tree.short

# rm family column from df
df<-subset(df, select = c( -Family))
dim(df) # 46 families
# make matrix

m <- df
summary(m)
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")
svg(file="heatmap_phyl_active.svg",width = 6, height=6)
windows(6,6)

p<- ggtree(tree.short, branch.length = .001) + 
 geom_tiplab(size=4, align=TRUE, linesize=0, offset = 0) + 
 theme_tree2()+
 xlim_tree(1) 
gheatmap(p, m, 
         colnames=FALSE,
         legend_title="active taxa", offset = 1.5, font.size = 10,) +
  scale_x_ggtree() + 
  scale_fill_gradient(low= "#fffaa2", high =  "#bb0000", aesthetics = "fill", na.value = "white",
                      name="Abundance in Active")+
  ggtitle("Families in Active Community")
dev.off()


##viable###
#select viable cell communtiy, ASVs with at least in 3 samples min reads is 50
df<-subset_samples(ps, Fraction=="Viable_Cell")
df<-ps_prune(df, min.samples = 3, min.reads = 50)
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
#remove "other" ASVS where rare taxa went
df<-filter(df, !grepl("Others", row.names(df)))
# remove "other" row where rare taxa went
taxon<-taxon[c(which(row.names(taxon)!="Others")),]


####### aggregate to the family level
df$otu<-row.names(df)
taxon$otu <- row.names(taxon)
df<-left_join(df, taxon)

#### change some the names so they match the tree
df$Family[grepl(" f__Xantho", df$Family)] <- " f__Xanthobacteraceae"  
df$Family[grepl(" f__Rhizo", df$Family)] <- " f__Rhizobiaceae"  
df$Family[grepl(" f__Pyrino", df$Family)] <- " f__Pyrinomonadaceae"
df$Family[grepl(" f__Burkhold", df$Family)] <- " f__Burkholderiaceae"
df$Family[grepl(" f__Solirub", df$Family)]  <-" f__Solirubrobacteraceae"
df$Family[grepl(" f__Chitino", df$Family)]  <-" f__Chitinophagaceae"
df$Family[grepl(" f__Rhodano", df$Family)] <- " f__Rhodanobacteraceae"
df$Family[grepl(" f__Blastocatellaceae", df$Family)] <- " f__Blastocatellaceae"

#rm other columns
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

colnames(df)
head(df)
#aggregate
df<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,
                    C1R.SYBR_S16,  C2E.SYBR_S22, C2N.SYBR_S15,  C2R.SYBR_S17,
                    C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25, C7R.SYBR_S19)
                     ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) #103
head(df)

#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.SYB"))) %>%   glimpse()
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)

#rm unknown families
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")
df<-filter(df, !grepl("UBA", df$otu))
df<-filter(df, !grepl("SG", df$otu))
df<-filter(df, !grepl("SCN", df$otu))
df$otu

# rename column
df<-mutate(df, Family=otu) %>% select(., -otu)

#removed 'f__'
df$Family<-sub(" f__", "", df$Family)
row.names(df) <- df$Family


#importtree
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/trees/GTDB")
tree = read.tree("family.nwk")
tree
# 1057 tips
#length(tree$tip.label) # look at the tip labels 
# modify tip labels
#tree$tip.label <- paste0(" f__", tree$tip.label)
length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
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
dim(df)
tree.short

# rm family column from df
df<-subset(df, select = c( -Family))
dim(df) # 49 families
# make matrix

m <- df
summary(m)
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m


setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures")
svg(file="heatmap_phy_viable.svg",width = 6, height=6)
windows(6,6)

p<- ggtree(tree.short, branch.length = .001) + 
  geom_tiplab(size=4, align=TRUE, linesize=0, offset = 0) + 
  theme_tree2()+
  xlim_tree(1) 
gheatmap(p, m, 
         colnames=FALSE,
         legend_title="Viable taxa", offset = 1.5, font.size = 10) +
  scale_x_ggtree() + 
  scale_fill_gradient(low= "#96bceb", high =  "#033d85", aesthetics = "fill", na.value = "white",
                      name="Abundance in Viable")+
  ggtitle("Families in Viable Community")
dev.off()




###BINOMIAL model####
# rrarefied data to be able to compare samples 
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/")
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

## Determine minimum available reads per sample ##
min.s<-min(rowSums(asvs.t))
min.s
### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
asvs.r<-rrarefy(asvs.t, min.s)
dim(asvs.t)
dim(asvs.r)

###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Viable_Cell', 'POS'='Active_Cell', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Fraction, metadat$Compartment))

##---make phyloseq object with rarefied data -------#
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

# taxa that are in the plant and are unassigned
#afb244d96b70a4c948b205f7f7eea5c5 259366
#9bf55fb48ef29780111f9e54dd204793  33352
df<-subset_samples(ps, Compartment=="Roots" | Compartment=="Nodule" )
df<-prune_taxa(taxa_sums(df) > 0, df)
remove<-subset_taxa(df, Domain=="Unassigned" |  Phyla=="" | Phyla==" p__"  ) 
remove
# 53 taxa
unique(taxon$Domain)
########## remove these guys
badtaxa<-taxa_names(remove)
alltaxa<-taxa_names(ps)
mytaxa <- alltaxa[!(alltaxa %in% badtaxa)]
ps<-prune_taxa(mytaxa, ps )
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
ps
# 12652 taxa

###make data frame##
#presence in plant ~ abundance in rhizosphere
# don't include taxa that are super rare less than 10
# in 3 samples

ps1<-subset_samples(ps, Compartment != "ctl"& Fraction != "Total_DNA")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 50)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1

# 927 asvs
taxon<- as.data.frame(tax_table(ps1))
df<-as.data.frame(otu_table(ps1))
#remove "other" column where rare taxa went
df<-select(df, -Others)
# remove "other" row where rare taxa went
taxon<-taxon[c(which(row.names(taxon)!="Others")),]

#can't use C67 (rep4) because it doesn't have active in the rhizosphere.

#### seperate out each rep
#               abundance_active  abundance_total present_in_plant
#otu #1 rep 1
#otu #1 rep 2
n<-row.names.data.frame(df)
n[grepl("C10" , n)]="5"
n[grepl("C2" , n)]="2"
n[grepl("C5" , n)]="3"
n[grepl("C7" , n)]="4"
n[grepl("C1R" , n)]="1"
n[grepl("C1N" , n)]="1"
n[grepl("C1E" , n)]="1"

df$rep <- n

#rep 1
df1<-filter(df, rep=="1")
dim(df1)
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
head(df1)
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
head(df1)
row.names(df1)<-paste0(row.names(df1), "Rep1")
colnames(df1) <- c("Active", "Total", "inplant")
rep1<-df1
head(rep1)

##rep 2
df1<-filter(df, rep=="2")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep2")
colnames(df1) <- c("Active", "Total", "inplant")
df1
rep2<-df1

##rep 3
df1<-filter(df, rep=="3")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep3")
colnames(df1) <- c("Active", "Total", "inplant")
df1
rep3<-df1

##rep 4 only has sybr
df1<-filter(df, rep=="4")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep4")
colnames(df1) <- c("Total", "inplant")
df1
rep4<-df1
dim(rep4)



##rep 5
df1<-filter(df, rep=="5")
#find out if somehting is in plant
n<-row.names.data.frame(df1)
n[grepl("N." , n)]="nodule"
n[grepl("E." , n)]="roots"
#remove rep
df1<-select(df1, -rep)
# sum by group and make vector
df1<-rowsum(df1, n)
df1<-as.data.frame(t(df1))
# insert column
inplant<-rep(NA, length(df1$nodule))
inplant[df1$nodule>50 | df1$roots>50 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep5")
#change column names
colnames(df1) <- c("Active", "Total", "inplant")
rep5<-df1
dim(rep5)
rep5
### put them together
df1<-rbind(rep1, rep2)  
df1<-rbind(df1, rep3)
df1<-rbind(df1, rep5)

df1$inplant <- as.numeric(df1$inplant)
colnames(df1)
dim(df1)
head(df1)


# Models:
fit <- glm(df1$inplant ~ df1$Total, family = binomial)
summary(fit)

fit <- glm(df1$inplant ~ df1$Total+ df1$Active + df1$Active*df1$Total, family = binomial)
summary(fit)

fit <- glm(df1$inplant ~ df1$Active , family = binomial)
summary(fit)

#plots :
ggplot(df1)+
  geom_point(aes(Total, inplant))+
  stat_smooth(aes(Total, inplant), method="glm", color="#045bc2", se=FALSE, 
                method.args = list(family=binomial))+
  xlim(0,3800)+
  theme_bw()


ggplot(df1)+
  geom_point(aes(Active, inplant))+
    stat_smooth(aes(Active, inplant), method="glm", color="#045bc2", se=FALSE, 
              method.args = list(family=binomial))+
  xlim(0,3800)+
  theme_bw()

##########TOP TAXA effecting the model #######
df1<-df1[order(-df1$Active),]
df2<-filter(df1, inplant==1)


# make rep column and average
head(df2)
otu<-row.names(df2)
rep<-str_sub(otu, -4, -1)

otu<-gsub("Rep1", "", otu)
otu<-gsub("Rep2", "", otu)
otu<-gsub("Rep3", "", otu)
otu<-gsub("Rep4", "", otu)
otu<-gsub("Rep5", "", otu)

df2$otu <- otu
df2$rep <- rep

# summarise by otu
df2<-df2 %>% group_by(otu) %>% summarise(mean_active = mean(Active), sd_active=sd(Active), mean_viable=mean(Total), sd_viable=sd(Total))

df2 <- filter(df2, mean_active>0)
df2

taxon$otu <-row.names(taxon) 
top_taxa<-left_join(df2, taxon)
head(top_taxa)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/tables")
write.csv(top_taxa, "activetaxa_inplant.csv")

