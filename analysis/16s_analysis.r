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

#basic
library(tidyverse)
library(vegan)
library(reshape2)
library(scales)
library(data.table)
#library(readxl)

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
#library(NatParksPalettes)
#library(viridis)
#library(hrbrthemes)
#library(ggpattern)
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


#-----load functions-------------#
# functions
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

#library(NatParksPalettes)
#names(NatParksPalettes)

#natparks.pals(name="Arches",n=8,type="discrete")
#natparks.pals(name="Arches2",n=4,type="discrete")
#natparks.pals(name="Banff",n=7,type="discrete")
#natparks.pals(name="Cuyahoga",n=6,type="discrete")


##########get colors#################
#from national park pal
mycols6<-c("#8fcafd",  "#4499f5" ,  "#0c62af" ,   "#f0ac7d", "#cd622e" , "#993203")

mycols8 <- c( "#8fcafd",  "#4499f5" ,  "#0c62af" ,   "#f0ac7d", "#cd622e" , "#993203", "#b46db3", "#3a1f46")
square <- 22
diamond <- 23
triangle <- 24
circle <- 21

mycols18<- c( "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A" ,
              "#FFFF99", "#B15928","#9ad6ce", "#1a635a", "#edceeb", "#eb05db", "#969696", "#232423")

phycols7<- c("#E31A1C", "#FB9A99","#FF7F00" , "#FDBF6F",  "#1F78B4","#A6CEE3", "#1a635a")

#####Import data ##########
## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/")
### Import Data silva###
#taxon <- read.table("asv_level_output/silva/taxonomy.txt", sep="\t", header=T, row.names=1)
#asvs.raw <- read.table("asv_level_output/silva/feature-table.tsv", sep="\t", header=T, row.names = 1 )
#metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)
#asvs.seq  = system.file("extdata", "asv_level_output/silva/dna-sequences.fasta", package="phyloseq")

### Import Data Greengenes###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)
#asvs.seq  = system.file("extdata", "asv_level_output/greengenes/dna-sequences.fasta", package="phyloseq")


## Transpose ASVS table ##
asvs.t <- t(asvs.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs.t<-asvs.t[order(row.names(asvs.t)),]
# make it percent
#asvs.perc<-((asvs.t/rowSums(asvs.t))*100)

#reads<-rowSums(asvs.t) 
#samples <- rownames(asvs.t)

#reads<-as.data.frame(cbind(samples, reads))
#row.names(reads) <-NULL
#reads
###--- recode metadata----- #
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

##------make phyloseq object with percent data -------#
asvs.phyloseq<- (asvs.t)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)


#import it phyloseq
Workshop_OTU <- otu_table(as.matrix(asvs.phyloseq), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833


ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat )
ps
# 15027 taxa silva
# 12855 taxa greengenes2

# who are the Otus that there taxa is unassigned?? 

df<-subset_taxa(ps, Phyla=="" | Phyla == " p__")
df<-as.data.frame(t(otu_table(df)))
df<-select(df, contains ("E"))
df1<-as.data.frame(rowSums(df))
colnames(df1) <- "sum"
df1<-filter(df1, sum!=0)
head(df1)
259366/sum(df1)
#sequece of 2 most abundant otus
# GACGGGGGGGGCAAGTGTTCTTCGGAATGACTGGGCGTAAAGGGCACGTAGGCGGTGAATCGGGTTGAAAGTTAAAGTCGCCAAAAACTGGTGGAATGCTCTCGAAACCAATTCACTTGAGTGAGACAGAGGAGAGTGGAATTTCGTGTGTAGGGGTGAAATCCGCAGATCTACGAAGGAACGCCAAAAGCGAAGGCAGCTCTCTGGGTCCCTACCGACGCTGGAGTGCGAAAGCATGGGGAGCGAACGGGA
# GACGGGGGGGGCAAGTGTTCTTCAGAATGACTGGGCGTAAAGGGCACGTAGGCGGTGAATCGGGTTGAAAGTTAAAGTCGCCAAAAACTGGTGGAATGCTCTCGAAACCAATTCACTTGAGTGAGACAGAGGAGAGTGGAATTTCGTGTGTAGGGGTGAAATCCGCAGATCTACGAAGGAACGCCAAAAGCGAAGGCAGCTCTCTGGGTCCCTACCGACGCTGGAGTGCGAAAGCATGGGGAGCGAACGGGA
# remove 
#afb244d96b70a4c948b205f7f7eea5c5 259366
#9bf55fb48ef29780111f9e54dd204793  33352

##------ percent abundance figure

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

#make rownames null
# summarize by phyla
df1<-aggregate(cbind(BEADS , Bulksoil.DNA.1  , Endo.BONCAT.1 , Nodule.BONCAT.1 , Nodule.Totalcells.1 , Rhizo.DNA.1,  Rhizo.BONCAT.1 , Rhizo.Totalcells.1 ,Endo.BONCAT.2 ,  Endo.Totalcells.2, 
Nodule.BONCAT.2 ,  Nodule.Totalcells.2  , Rhizo.DNA.2,   Rhizo.BONCAT.2,   Rhizo.Totalcells.2,  Bulksoil.DNA.3,    Endo.BONCAT.3,   Endo.Totalcells.3 , Nodule.BONCAT.3,   Nodule.Totalcells.3 ,
Rhizo.DNA.3 ,  Rhizo.BONCAT.3  , Rhizo.Totalcells.3 , Bulksoil.DNA.4 ,    Endo.BONCAT.4,    Endo.Totalcells.4,  Nodule.BONCAT.4 ,  Rhizo.DNA.4,   Rhizo.BONCAT.4,   Rhizo.Totalcells.4 ,
 BulkSoil.DNA.5,   Endo.BONCAT.5,   Endo.Totalcells.5,  Nodule.BONCAT.5 ,  Nodule.Totalcells.5 , Rhizo.DNA.5,   Rhizo.Totalcells.5 , CTL  ,     Bulksoil.DNA.6  ,  Bulksoil.DNA.7  ,  
 Bulksoil.DNA.8,     Bulksoil.DNA.9) ~ Phyla, data = df1, FUN = sum, na.rm = TRUE)

# call empty phyla unassigned
df1$Phyla[c(1,2)] <- "Unassigned"
head(df1)

# gather by sample
df1<-gather(df1, "sample", value, 2:43 )
head(df1)
#remove zeros
df1<-df1[df1$value!=0,]
head(df1)
#df1<-df # place holder
# order columns

# make really low abundance taxa other
df1$Phyla[df1$value<1] <- "other"
df1<-aggregate(cbind(value) ~ sample+Phyla, data = df1, FUN = sum, na.rm =TRUE)
head(df1)
df1<-df1[order(df1$sample),]

windows(12,10)
df1%>% 
  ggplot(aes(fill=Phyla, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values=mycols18) +
  #scale_fill_viridis(discrete = TRUE) +
  ggtitle("Top phyla") +
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

dev.off()

###### remove unassigned taxa  ########
# taxa that are in the plant and are unassigned
# remove 
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
# 12802 taxa



##### percent abundance figure No. 2######
###### remake figure
# grab data
taxon<- as.data.frame(tax_table(ps))
df<- as.data.frame(otu_table(ps))
# make it percent
df<-(df/rowSums(df))*100
df<-as.data.frame(t(df))
df<-cbind(df, taxon)

head(df)
# summarize by phyla
#make rownames null
n<-c("BEADS" ,    "Bulksoil.DNA.1"  , "Endo.BONCAT.1" , "Nodule.BONCAT.1" , "Nodule.Totalcells.1" ,"Rhizo.DNA.1",  "Rhizo.BONCAT.1" , "Rhizo.Totalcells.1" ,"Endo.BONCAT.2" ,  "Endo.Totalcells.2", 
     "Nodule.BONCAT.2" ,  "Nodule.Totalcells.2"  ,"Rhizo.DNA.2" ,   "Rhizo.BONCAT.2",   "Rhizo.Totalcells.2",  "Bulksoil.DNA.3",    "Endo.BONCAT.3",   "Endo.Totalcells.3" , "Nodule.BONCAT.3",   "Nodule.Totalcells.3" ,
     "Rhizo.DNA.3" ,  "Rhizo.BONCAT.3"  , "Rhizo.Totalcells.3" , "Bulksoil.DNA.4" ,    "Endo.BONCAT.4",    "Endo.Totalcells.4",  "Nodule.BONCAT.4" ,  "Rhizo.DNA.4",   "Rhizo.BONCAT.4",   "Rhizo.Totalcells.4" ,
     "BulkSoil.DNA.5",   "Endo.BONCAT.5",   "Endo.Totalcells.5",  "Nodule.BONCAT.5" ,  "Nodule.Totalcells.5" , "Rhizo.DNA.5",   "Rhizo.Totalcells.5" , "CTL"  ,     "Bulksoil.DNA.6"  ,  "Bulksoil.DNA.7"  ,  
     "Bulksoil.DNA.8",     "Bulksoil.DNA.9"  ,   "Domain"  ,      "Phyla"       ,  "Class"  ,       "Order"    ,     "Family"  ,      "Genus"    ,     "Species"   )
df<-df
colnames(df)<-n
row.names(df) <- NULL

# summarize by phyla
df<-aggregate(cbind(BEADS , Bulksoil.DNA.1  , Endo.BONCAT.1 , Nodule.BONCAT.1 , Nodule.Totalcells.1 , Rhizo.DNA.1,  Rhizo.BONCAT.1 , Rhizo.Totalcells.1 ,Endo.BONCAT.2 ,  Endo.Totalcells.2, 
                     Nodule.BONCAT.2 ,  Nodule.Totalcells.2  , Rhizo.DNA.2,   Rhizo.BONCAT.2,   Rhizo.Totalcells.2,  Bulksoil.DNA.3,    Endo.BONCAT.3,   Endo.Totalcells.3 , Nodule.BONCAT.3,   Nodule.Totalcells.3 ,
                     Rhizo.DNA.3 ,  Rhizo.BONCAT.3  , Rhizo.Totalcells.3 , Bulksoil.DNA.4 ,    Endo.BONCAT.4,    Endo.Totalcells.4,  Nodule.BONCAT.4 ,  Rhizo.DNA.4,   Rhizo.BONCAT.4,   Rhizo.Totalcells.4 ,
                     BulkSoil.DNA.5,   Endo.BONCAT.5,   Endo.Totalcells.5,  Nodule.BONCAT.5 ,  Nodule.Totalcells.5 , Rhizo.DNA.5,   Rhizo.Totalcells.5 , CTL  ,     Bulksoil.DNA.6  ,  Bulksoil.DNA.7  ,  
                     Bulksoil.DNA.8,     Bulksoil.DNA.9) ~ Phyla, data = df, FUN = sum, na.rm = TRUE)




tail(df)

#gather
df<-gather(df, "sample", value, 2:43 )
df
#remove zeros
# anything that is less than 1% equals other
df<-df[df$value!=0,]

#df1<-df # place holder
df$Phyla[df$value<1] <- "other"
unique(df$Phyla)


mycols17<- c( "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A" ,
              "#FFFF99", "#B15928","#9ad6ce", "#1a635a", "#edceeb", "#eb05db", "#969696")

df<-df[order(df$sample),]



windows(10,10)
df%>% 
  ggplot(aes(fill=Phyla, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values = mycols18)+
  #scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  #scale_fill_viridis(discrete = TRUE) +
  ggtitle("Top phyla") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

# use this ps!
ps

# this but just rhizosphere

filter(df, sample!="BEADS")
df1<-df[grep("Rhizo" , df$sample),]
Total<-df1[grep("Total" , df1$sample),]
bcat<-df1[grep("BONCAT" , df1$sample),]
df<-rbind(Total, bcat)
df<-df[order(df$sample),]
df$Phyla[df$value<1] <- "other"
unique(df$Phyla)

windows(10,10)
df%>% 
  ggplot(aes(fill=Phyla, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values = mycols17)+
  #scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  #scale_fill_viridis(discrete = TRUE) +
  ggtitle("Top phyla") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
dev.off()

# mak ea plot liek ryans



########barplot of endophyte taxa

endo<-subset_samples(ps, Compartment == "Roots")
endo<-prune_taxa(taxa_sums(endo) > 0, endo)
endo
# 648 taxa

# grab data
taxon<- as.data.frame(tax_table(endo))
df<- as.data.frame(otu_table(endo))
# make it percent
df<-(df/rowSums(df))*100
df<-as.data.frame(t(df))
df<-cbind(df, taxon)

# summarize by phyla
df<-aggregate(cbind(C10E.POS_S60  , C1E.POS_S31,   C1E.SYBR_S21,  C2E.POS_S32,   C2E.SYBR_S22,  C5E.POS_S33,
                    C5E.SYBR_S23, C7E.POS_S34,   C7E.SYBR_S24) ~ Genus, data = df, FUN = sum, na.rm = TRUE)
head(df)
df$Genus[c(1,2)] <- "other"

#gather
df<-gather(df, "sample", value, 2:10 )
#remove zeros
df<-df[df$value!=0,]


windows(10,10)
df%>% 
  ggplot(aes(fill=Genus, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  #scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Genus") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.position="none")
dev.off()

df


##old import Asvs data ##
## Set the working directory; ###
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")

### Import Data ###
#taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
#asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
#metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)

## Transpose ASVS table ##
#asvs.t <- t(asvs.raw)
## order metadata
#metadat<-metadat[order(metadat$SampleID),]
## order asvs table
#asvs.t<-asvs.t[order(row.names(asvs.t)),]

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
#metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
#metadat<-metadat[, c(1,3:6)]
#metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
#metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

##------make phyloseq object with rarefied data -------#

#row.names(df) <- df$column1


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

#test it worked
#sample_names(ps)
#print(ps)
# 12855 taxa

# remove chloroplast DNA
#ps<-subset_taxa(ps, Class!=" Chloroplast")
#ps<-subset_taxa(ps, Genus!=" Mitochondria")
#ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that aren; in any samples
#ps<-prune_taxa(taxa_sums(ps) > 0, ps)
#any(taxa_sums(ps) == 0)
#ps
# 12855 taxa

#asvs.clean<-as.data.frame(t(as.data.frame(otu_table(ps))))
#taxon<-as.data.frame(tax_table(ps))



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


# 3173 taxa
# output df 
#asvs.clean<-as.data.frame(t(as.data.frame(otu_table(ps))))
#taxon<-as.data.frame(tax_table(ps))



#####DIVERSITY  figure#####
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
#how many ASVS in each fraction?
rich %>% group_by(rich$compartment_BCAT) %>% summarise(mean(Observed), sd(Observed))
rich %>% group_by(rich$Compartment) %>% summarise(mean(Observed), sd(Observed))


# total + active number otus
svg(file="figures/16s/observed_divserity.svg",width = 10, height=6 )
windows(width = 6, height=3.5)
rich %>%
  ggplot(aes(x=Compartment, y=Observed,  col= Fraction))+
  geom_boxplot() +
  scale_color_manual(values=mycols8[c(7,1,5)]) +
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")+
  facet_wrap(~Fraction, scales = "free_x")+
  scale_x_discrete(drop = TRUE) +
  ylab("N asvs")+
  xlab("Compartment")
dev.off()

# total + active number otus
svg(file="figures/16s/shannon_divserity.svg",width = 8, height=5 )
windows(width = 8, height=5)
rich %>%
  ggplot(aes(x=Compartment, y=Shannon,  col= Fraction))+
  geom_boxplot() +
  scale_color_manual(values=mycols8[c(7,1,5)]) +
  geom_jitter(width = .1, size=1 )+
  theme_classic(base_size = 18)+
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")+
  facet_wrap(~Fraction, scales = "free_x")+
  scale_x_discrete(drop = TRUE) +
  ylab("Shannon")+
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






######DIVERSITY Stats######
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

root<-subset_samples(ps, Compartment=="Roots")
root<-prune_taxa(taxa_sums(root) > 0, root)
tax_table(root)
root



## fitler for taxa that are rare
root_total<-prune_taxa(taxa_sums(root_total) < 10, root_total)
any(taxa_sums(root_total) == 0)
root_total
# 231 taxa

x<-otu_table(root_total)
x
#calc richness
rich<-estimate_richness(root, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

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
# are taxa in the root present in the bulk soil ###
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


#####ESTELLE fig. BARPLOT of active verse total abundance for top otus######

#BAR CHART with TOP 50 OTU##
# filter for total and active
ps
sample_data(ps)
df<-subset_samples(ps, Compartment!="ctl"& Compartment!="Bulksoil" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df<-as.data.frame(t(otu_table(df)))
#dim(df)
#df
# 5396 taxa
df<-cbind(taxon,df)

#get means
df$nodule.total.mean <-   rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
df$nodule.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root.total.mean <-   rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
df$root.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()

colnames(df)

# select top 50 otus rhizosphere
otu<-df%>% select(-contains("POS")) %>% select(-contains("SYBR")) 
head(otu)
otu<-otu[order(otu$rhizo.total.mean, decreasing = TRUE),]
top<-otu[1:50,]
top$otu<-row.names(top)
top

# what percent of the whole community are the top asvs
all<-colSums(otu[,12:13])
t<-colSums(top[,12:13])
t/all
# top 50 otus is 34% of the total population;
write.csv(top, file="top50asvs.csv")

# lil labels
df <- data.frame(x     = (top$otu ),
                 y     = (top$rhizo.bcat.mean + 300) ,
                 label = c("1", "", "2*", "3", rep("",8 ), "4", "",
                           "5", "", "6", rep("",16), "7", rep("", 8), "8*", rep("", 7) ))
df
top
## 1  Rhizobium
## 2** Sphingomonadaceae
## 3  Pseudomonas_E_647464                                              
## 4 Acidimicrobiales
## 5 Chthoniobacteraceae
## 6 Actinobacteria family Micrococcaceae
## 7 Rhizobiaceae
## 8** DA Sphingomonadaceae


#edit labels
top$Phyla<-sub("p__", "", top$Phyla)
top<-top[order(top$rhizo.total.mean, decreasing = TRUE),]
top$otu1<-c(1:50)
top
#plot
windows(6,5)
ggplot(top)+
  geom_bar(aes(x=reorder(otu1, -rhizo.total.mean) , y=rhizo.bcat.mean, fill= Phyla), stat="identity", position="identity") +
  geom_line(aes(x=as.numeric(reorder(otu1, -rhizo.total.mean)) , y=rhizo.total.mean), linewidth=1) +
  xlab("top 50 otus")+
  ylab("Abundance")+
  scale_fill_manual(values=phycols7)+
  theme_classic(base_size = 16)+
  theme(axis.text.x=element_blank(), legend.position = c(.65, .7))+
  geom_text(data = df, aes(x, y, label = label), size=5)+
  guides(fill=guide_legend(title="Abundance Active Cells"))

#legend.position = c(.65, .7)
# in plant is red.    
ggplot(bcat, aes(x=reorder(otu, +order) , y=log_abundance, fill = location, pattern = location)) +
    geom_bar_pattern(stat="identity", position = "dodge",
                     color = "black", 
                     pattern_fill = "black",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values=   c("black","grey", "red", "red"))+
  scale_pattern_manual(values = c(total_rhizo = "none",nodule = "none", rhizo = "none", root = "stripe")) +
  labs(x = "otu", y = "abundance", pattern = "location") +
  theme_bw()
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))  
    
#####ESTELLE fig. BARPLOT but abundance in endophyte is there too. ######
    
    #BAR CHART with TOP 50 OTU##
    # filter for total and active
    ps
    sample_data(ps)
    df<-subset_samples(ps, Compartment!="ctl"& Compartment!="Bulksoil" & Fraction!="Total_DNA")
    df<-prune_taxa(taxa_sums(df) > 0, df)
    taxon<-tax_table(df)
    df<-as.data.frame(t(otu_table(df)))
    #dim(df)
    #df
    # 5396 taxa
    df<-cbind(taxon,df)
    
    head(df)
    #nodule
    df$nodule.total.mean <-   rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
    #df$nodule.total.log <- rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>% log10() %>%  glimpse()
    df$nodule.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
    #df$nodule.xbcat.log <- rowMeans(df %>% dplyr::select(contains("N.POS"))) %>% log10() %>%  glimpse()
    # roots
    df$root.total.mean <-   rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
    #df$root.total.log <- rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% log10() %>%  glimpse()
    df$root.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
    #df$root.xbcat.log <- rowMeans(df %>%  dplyr::select(contains("E.POS"))) %>% log10() %>%  glimpse()
    
    #rhizo
    df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
    # add 1 to everything so there are no negative number
    #df$rhizo.total.log<-(df$rhizo.total.mean +1 ) %>% log10()
    df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
    # add 1 to everything so there are no negative number
    #df$rhizo.bcat.log<-(df$rhizo.bcat.mean +1 ) %>% log10()
    
    #df[df== "-Inf"] <- 0
    colnames(df)
    
    # select taxa that are present in the plant 
    otu<-df%>% select(-contains("POS")) %>% select(-contains("SYBR")) 
    otu<-otu[otu$root.xbcat.mean>0 &  otu$root.total.mean >0 & otu$nodule.xbcat.mean>0 & otu$nodule.total.mean>0  & otu$rhizo.total.mean>0 & otu$rhizo.bcat.mean>0,]
    dim(otu)    # 9 taxa present everywhere
    # sort by abundance
    otu<-otu[order(otu$root.total.mean, decreasing = TRUE),]
    otu
    
         write.csv(otu, file="topgeneralistasvs.csv")
    # lil labels
    df <- data.frame(x     = top$otu,
                     y     = (top$rhizo.bcat.mean + 70) ,
                     label = c("*1", "*2", "*3", "*4", "*5", "*6", rep("", 44)))
    df
    
    #plot
    windows(6,5)
    ggplot(top)+
      geom_bar(aes(x=reorder(otu, -rhizo.total.mean) , y=rhizo.bcat.mean, fill= Phyla), stat="identity", position="identity") +
      geom_line(aes(x=as.numeric(reorder(otu, -rhizo.total.mean)) , y=rhizo.total.mean)) +
      xlab("top 50 otus")+
      ylab("Abundance")+
      scale_fill_manual(values=phycols7)+
      theme_classic(base_size = 16)+
      theme(axis.text.x=element_blank(), legend.position = c(.65, .7) )+
      geom_text(data = df, aes(x, y, label = label))+
      guides(fill=guide_legend(title="Abundance Active Cells"))
    
    
    

# in plant is red.    
ggplot(bcat, aes(x=reorder(otu, +order) , y=log_abundance, fill = location, pattern = location)) +
    geom_bar_pattern(stat="identity", position = "dodge",
                     color = "black", 
                     pattern_fill = "black",
                     pattern_angle = 45,
                     pattern_density = 0.1,
                     pattern_spacing = 0.025,
                     pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values=   c("black","grey", "red", "red"))+
  scale_pattern_manual(values = c(total_rhizo = "none",nodule = "none", rhizo = "none", root = "stripe")) +
  labs(x = "otu", y = "abundance", pattern = "location") +
  theme_bw()
    guides(pattern = guide_legend(override.aes = list(fill = "white")),
           fill = guide_legend(override.aes = list(pattern = "none")))    
    

#####RYAN plot shifts in from total to active phyla level#####
# filter for total and active
    # filter for total and active
    ps
    sample_data(ps)
    df<-subset_samples(ps, Compartment!="ctl"& Compartment=="Rhizosphere" & Fraction!="Total_DNA")
    df<-prune_taxa(taxa_sums(df) > 0, df)
    taxon<-tax_table(df)
    df<-as.data.frame(otu_table(df))
    dim(df)
    #df
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
#df[ , order(names(df))]
top<-top%>%
  pivot_longer(2:10, values_to = "abundance", names_to= "rep_fraction" )
top<-top%>%mutate(fraction = str_split_i(top$rep_fraction, "_", 1)) %>%
  mutate(Phyla= str_split_i(top$Phyla, "_", 3))

head(top)
# ploy
top$fraction<-factor(top$fraction, levels = c("Total", "BONCAT"))

windows(7,3)
ggplot(top) +
  geom_boxplot(aes(x= reorder(fraction, -sum), y= abundance, fill=Phyla), outlier.shape = NA, alpha=.5, show.legend = FALSE)+
  geom_jitter(aes(x= reorder(fraction, -sum), y= abundance, fill=Phyla), position = position_jitter(width = .2), show.legend = FALSE)+
  facet_wrap(~reorder(Phyla, -sum), nrow=1)+
  theme_classic(base_size = 12)+
  scale_fill_manual(values=phycols7)+
  labs(y="proportion of reads in sample")+
  theme(axis.text.x=element_text(angle=45, hjust=0.9), axis.title.x = element_blank())




######################### ANCOM ###############
sample_data(ps)
ps1<-subset_samples(ps, Compartment=="Rhizosphere" & BONCAT!="DNA")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 100)
sample_data(ps1)
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
#
# add otu column
taxon<-as.data.frame(tax_table(ps1))
taxon$otu<-row.names(taxon)
tax_table(ps1) <- tax_table(as.matrix(taxon))
tax_table(ps1)
ps1
# 5250 taxa if you want to rm taxa that are super rare
sample_data(ps1)

ps2 = mia::makeTreeSummarizedExperimentFromPhyloseq(ps1)

ps2

out1 = ancombc(data = ps2, assay_name = "counts", 
               tax_level = "otu", phyloseq = NULL, 
               formula = "Fraction", 
               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
               group = "Fraction", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
               n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
sample_data(ps1)
# ancombc also supports importing data in phyloseq format
# tse_alt = agglomerateByRank(tse, "Family")
# pseq = makePhyloseqFromTreeSummarizedE


#log fold change
tab_lfc = res$lfc
head(tab_lfc)
dim(tab_lfc)
col_name = c("otu", "LFC_Intercept", "LFC_Fraction_Total_cells_Active")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
head(tab_lfc)
#write_delim(as.data.frame(tab_lfc), file = "Log_fold_change_ancom.txt", delim = " ")

# stadnard error
tab_se = res$se
x<-tab_se$FractionTotal_Cells
x<-as.numeric(x)
names(x)<-NULL
head(x)

otu<-tab_se$taxon
head(otu)
tab_se <- data.frame (otu = otu,
                  SE_Fraction_Total_cells_Active = x
)
head(tab_se)
#write_delim(as.data.frame(tab_se), file = "SE_ancom.txt", delim = " ")

#test statistcs W maybe it's willcoxin?
tab_w = res$W
col_name = c("otu", "W_Intercept", "W_Fraction_Total_cells_Active")
colnames(tab_w) = col_name
head(tab_w)
#write_delim(as.data.frame(tab_se), file = "SE_ancom.txt", delim = " ")

# "P-values from the Primary Result
tab_p = res$p_val
x<-tab_p$FractionTotal_Cells
x<-as.numeric(x)
names(x)<-NULL
otu<-tab_p$taxon
head(otu)
tab_p <- data.frame (otu = otu,
                     pval_Fraction_Total_cells_Active= x
)
head(tab_p)
#write_delim(as.data.frame(tab_), file = "pval_ancom.txt", delim = " ")

#Adjusted p-values from the Primary Result"
tab_q = res$q
x<-tab_q$FractionTotal_Cells
x<-as.numeric(x)
names(x)<-NULL
otu<-tab_q$taxon
tab_q <- data.frame (otu = otu,
                     adj_pval_Fraction_Total_cells_Active= x
)
head(tab_q)
#write_delim(as.data.frame(tab_se), file = "adjpval_ancom.txt", delim = " ")

# yes or no is a taxa differentially abundant
tab_diff = res$diff_abn
col_name = c("otu", "DA_Intercept", "DA_Fraction_Total_cells_Active")
colnames(tab_diff) = col_name
head(tab_diff)
#write_delim(as.data.frame(tab_se), file = "DA_ancom.txt", delim = " ")



#####Format df from ANCOM######
# through all togetha nd remove anything with DNA

tab <-tab_lfc %>%
    left_join(., tab_se) %>%
    left_join(., tab_w ) %>%
    left_join(., tab_p) %>%
    left_join(., tab_q) %>%
    left_join(., tab_diff)
write_csv(as.data.frame(tab), file = "table_ancom.csv", col_names = FALSE)

##add otu values to the output df
df<-as.data.frame(t(otu_table(ps1)))
head(df)
df$otu<-row.names(df)
df<-left_join(df, tab)
# taxon table
taxon<-as.data.frame(tax_table(ps1))
df<-left_join(taxon, df)

# filter for taxa that were not detected in the rhizosphere total 
df<-filter(df, df$C10R.POS_S30 > 0 | df$C1R.SYBR_S16 > 0 | df$C2R.SYBR_S17 > 0 | df$C5R.SYBR_S18 >0  )
dim(df) # 587 taxa

#make my own threshold too set threshold
#df2<-df1 %>% filter(pvalue_total_active<.01, abs(lfc_total_active)>1)
DA_taxa <- rep(NA, length(df$otu))
DA_taxa<-df$pval_Fraction_Total_cells_Active<.01 & abs(df$LFC_Fraction_Total_cells_Active)>1
DA_taxa
df$DA_taxa <- DA_taxa
#
length(which(df$DA_taxa=="TRUE"))
# 162 differential abundant taxa rhizosphere. 
length(which(df$DA_Fraction_Total_cells_Active=="TRUE"))
# 206 taxa
# I'm going to use both  my taxa threshold and the ancom onebecause it's a little more strict

# plotssss
# lfc verse the pvalue
ggplot(df, aes(x=LFC_Fraction_Total_cells_Active , y=pval_Fraction_Total_cells_Active)) + 
  geom_jitter()+ #(stat="identity", position="identity") +
  #scale_fill_manual(values=c("#ab8af2" ,"#4c4b4d"))+
  #xlab("differentially abundant otus")+
  theme_bw()

#plot
ggplot(df, aes(x=LFC_Fraction_Total_cells_Active , y=pval_Fraction_Total_cells_Active, col=DA_Fraction_Total_cells_Active)) + 
  geom_jitter()+ #(stat="identity", position="identity") +
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  theme_bw(base_size = 16)
#plot
ggplot(df, aes(x=LFC_Fraction_Total_cells_Active , y=-log(pval_Fraction_Total_cells_Active), col=DA_taxa)) + 
  geom_jitter()+ #(stat="identity", position="identity") +
  scale_color_manual(values=c("#999999", "#56B4E9"))+
  theme_bw(base_size = 16)

#LOOK AT DA taxa
# select taxa
DA<-df %>% filter(df$DA_taxa==TRUE) %>% filter(DA_Fraction_Total_cells_Active==TRUE)
dim(DA)
# 150 taxa match both
head(DA)
# add means
DA$rhizo.total.mean <-   rowMeans(DA %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
DA$rhizo.bcat.mean <-   rowMeans(DA %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
colnames(DA)

#rm some columsn to make things simpler
DA<-DA %>% select(-contains("R."))
dim(DA)
head(DA)

#save DA taxa file
write.csv(DA, file="DA_taxa.csv")

#order by abundance
DA<-DA[order(DA$rhizo.total.mean, decreasing = TRUE),]
head(DA)
top10<-DA[1:10,]


# multiply by -1 to make negative is lower in total
DA$LFC_Fraction_Total_cells_Active<-DA$LFC_Fraction_Total_cells_Active*-1

# sort for negative slope - more dormant taxa
dormant<-DA[DA$LFC_Fraction_Total_cells_Active < 0,]
dormant<-dormant[order(dormant$rhizo.total.mean, decreasing = TRUE),]
top10<-dormant[1:10,]
top10
write.csv(top10, file ="top10dormantaxa.csv")


# sort for more active taxa - postive slope
active<-DA[DA$LFC_Fraction_Total_cells_Active > 0,]

active<-active[order(active$rhizo.bcat.mean, decreasing = TRUE),]
top10_active<-active[1:10,]
top10_active

######### ANCOM - out of DA taxa how many are in the top 50 asvs ######
setwd("C:/Users/Jenn/The Pennsylvania State Universi#ty/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s")
top<-read_csv("top50asvs.csv")
top<-top[,2:5]
colnames(top) <- c("otu", "abudance", "label", "order")
head(top)

DA<-read_csv("DA_taxa.csv", row)
DA<-DA[,2:21]

both<-inner_join(DA, top)

###########phylo plot of differentially abundant taxa and there mean abundance#####

library(ggtree)
library(TDbook)
library(ape)
library(readxl)
library('TreeTools')

# grab tree
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree <- read.tree("trees/SEPPoutput/tree.nwk" )
tree$tip.label
taxon$Species

# make tree small to match differently abundant taxa.
df2
#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, df2$otu) #asvs we don't want
tree<-ape::drop.tip(tree, asvs_remove,  collapse.singles = TRUE,) # remove asvs we don't need
length(tree.short$tip.label)
tree
tree_seq_nwk
# grab correct order 
target<-tree.short$tip.label
df2<-df2[match(target, df2$otu),]
dim(df2)

plot(tree)
## load `tree_nwk`, `df_info`, `df_alleles`, and `df_bar_data` from 'TDbook'

colnames(df2)
## visualize the tree 

length(tree.short$tip.label)
p <- ggtree(tree.short) 
p
## attach the sampling information data set 
## and add symbols colored by location
#p <- p %<+% df2 + geom_tippoint(aes(color=Phyla))

## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
#grDevices::windows(14,10)
tree.short
head(df2)

p + geom_facet(panel = "Log 10 abundance in Active", data = df2, geom = geom_col, 
               aes(x=lfc_total_active, orientation = y, width = 1))
# geom_facet(panel = "Phyla", data = df2, geom = geom_text, size= 3 ,
#                       aes(x=0,label=ifelse(df2$first==1, unique(df2$Phyla), "")))+
#theme(legend.position="none")


dev.off()
# add node points
# geom_nodelab(node='internal', size=2, color="green", nudge_x = 1)+
#geom_cladelabel(node=600, label="Some random clade")
#?nodelabels
xc
# get hte first occurance
t.first <- match(unique(df2$Phyla), df2$Phyla)
t.first
df2[t.first,]
unique(df2$Phyla)

df2$first<- NULL
df2$first<-ifelse(match(row.names(df2), t.first), 1)

####################################3
plot(tree.short, show.tip.label = FALSE)+
  nodelabels(node=tree.short$edge[1172:1175,1])

n<-c(595, 602, 603, 610, 613, 614, 615, 616, 617, 674, 671, 750, 765, 766, 768, 844, 845, 846, 965, 1163, 1182 )

print(nodelabels(tree.short))
tree.short$edge[101:200,1]
#geom_cladelabel(node=17, label="Some random clade", 
#         
#color="red2", offset=.8, align=TRUE) + 

# nodelabels(text = tree.short$node.label,
#           frame = "n", cex=0.8, col= "blue")        

length(tree.short$edge[,1])

#####################are many taxa present in total cells not active ########

# how many otus are active
r<-df%>%select(contains("rhizo")) %>% filter(rhizo.bcat.sum>1)
dim(r)
# 2853 otus are active
hist(r$rhizo.bcat.sum, breaks = 50)
hist(log10(r$rhizo.bcat.sum), breaks = 50)
summary(r$rhizo.bcat.sum)
# most are pretty rare.
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     2.0    16.0    38.0   103.6    85.0 17518.0 

#how many are active but not in total cell.s
r<-df%>%select(contains("rhizo")) %>% filter(rhizo.total.sum<1)%>% filter(rhizo.bcat.sum>1)
dim(r)
# 1426

# how many taxa are active but not in total dna
r<-df%>%select(contains("rhizo")) %>% filter(rhizo.DNA.sum<1)%>% filter(rhizo.bcat.sum>1)
dim(r)
# 1706

### there are many taxa present in total dna not active
r<-df%>%select(contains("rhizo")) %>% filter(rhizo.DNA.sum>1)%>% filter(rhizo.bcat.sum<1)
dim(r)
# 2953  


r<-df%>%select(contains("rhizo")) %>% filter(rhizo.total.sum>1)%>% filter(rhizo.bcat.sum<1)
dim(r)
# 2399  
summary(r$rhizo.total.sum)
# typically rare taxa

#nodulesa
r<-df%>%select(contains("nodule")) %>% filter(nodule.xbcat.sum>0)
dim(r)

# 96 are active


r<-df%>%select(contains("nodule")) %>% filter(nodule.total.sum>0)
dim(r)
r
# 115 in total


r<-df%>%select(contains("nodule")) %>% filter(nodule.total.sum>0) %>% filter(nodule.xbcat.sum>0)
dim(r)
# 25 in both



#how many are active but not in total cell.s
r<-df%>%select(contains("nodule")) %>%filter(nodule.xbcat.sum>0) %>% filter(nodule.total.sum<1)
dim(r)
# 71

### there are many taxa present in total cells not active
r<-df%>%select(contains("nodule")) %>%filter(nodule.xbcat.sum<1) %>% filter(nodule.total.sum>1)
dim(r)
# 90  

# typically rare taxa


#######SCATTERPLOT active and total##############
sample_data(ps)
df<-subset_samples(ps, Compartment!="ctl" & Compartment=="Rhizosphere" & Fraction!="Total_DNA")
df<-prune_taxa(taxa_sums(df) > 0, df)
taxon<-tax_table(df)
df<-as.data.frame(t(otu_table(df)))
dim(df)
# 4854 taxa
#df<-cbind(taxon,df)

colnames(df)
n<-c( "BONCAT_1", "Total_1", "BONCAT_2" ,  "Total_2",  "BONCAT_3" , 
 "Total_3" , "BONCAT_4" ,  "Total_4",  "Total_5" )
colnames(df)<-n
#make a column for the value in active and a columsn for the value in total
#make a columns that stay the rep. 
df$otu <- row.names(df) 
total<-df %>% select(contains("Total"))
total$otu <- row.names(total) 
total<-total %>% pivot_longer(cols = 1:4, values_to = "total", names_to = "rep_total") %>% select(-Total_5)
#active
active<-df %>% select(contains("BONCAT"))
active$otu <- row.names(active) 
active<-active %>% pivot_longer(cols = 1:4, values_to = "active", names_to = "rep")
active[,2:3]
df<-cbind(total, active[,2:3])
df<-df%>% filter(total>10 & active>10 )

df<-df%>%
mutate( rep= str_split_i(df$rep, "_", 2))

#plot
windows(6,3)
filter(df) %>%
  ggplot(aes(x=log10(total) , y=log10(active), col=rep)) + 
  geom_jitter(stat="identity", position="identity") +
  #geom_abline(slope=1)+
  geom_smooth(method=lm, se=FALSE)+
  scale_color_natparks_d("KingsCanyon")+
  theme_bw(base_size = 12)+
  labs(x="Abundance in Total Viable Cell population",y="Abundance in Active population")
#lm
m1<-lm(log10(total)~log10(active), data=df)
summary(m1)
#R squared = .32

  

df%>%
  ggplot(aes(x=log10(nodule.xbcat.mean) , y=log10(nodule.total.mean))) + 
  geom_jitter(stat="identity", position="identity") +
  geom_abline(slope=1)+
  theme_bw()

 df%>%
   ggplot(aes(x=log10(root.xbcat.mean) , y=log10(root.total.mean))) + 
   geom_jitter(stat="identity", position="identity") +
   geom_abline(slope=1)+
   theme_bw()
 


#####VENN DIAGRAM  total#####
df<-subset_samples(ps, Fraction=="Total_Cells")
# at least in 2 samples min reads is 10
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df<-prune_taxa(taxa_sums(df) > 0, df)
df
dim(df)
# 1293 taxa


# don't include taxa that are super rare less then 10 reads
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))
#df<-df[-which(row.names(df)=="Others"),]

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
df<-as.data.frame(t(df))

df<-df[-which(row.names(df)=="Others"),]



dim(df)
habitat<-rep(NA, length(df$nodule))
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="rhizo specailist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="root specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="rhizo and nodule generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
unique(habitat)

### alternative labeling
habitat<-rep(NA, length(df$nodule))
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="generalist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
unique(habitat)


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
#x

################################## venn diagram #
#library(ggvenn)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s")
svg(file="OTU_level_total_venn.svg",width = 4, height=4 )
windows(4,4)
ggvenn(
  x, 
  fill_color = mycols6[c(1,2,3)],
  stroke_size = 2, set_name_size = 4, text_size = 3, digits = 1, fill_alpha=.6
  #auto_scale = TRUE
) +
  ggtitle("Total Cells")
dev.off()
#####VENN DIAGRAM active######
### don't include taxa that are only in 1 samples and less than 50 reads total
sample_data(ps)
# at least in 2 samples min reads is 10
df<-subset_samples(ps, Fraction=="BONCAT_Active")
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df<-prune_taxa(taxa_sums(df) > 0, df)
df
# 842 taxa
# rm others
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))
df
n<-row.names(df)
n[grepl("R" , n)]="rhizo"
n[grepl("E" , n)]="roots"
n[grepl("N" , n)]="nodule"
df
n
# summ by group
df<-rowsum(df, n)
df<-as.data.frame(t(df))
df

df<-df[-which(row.names(df)=="Others"),]

dim(df)
head(df)
habitat<-rep(NA, 805)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>1]="rhizo specailist"
habitat[df$nodule==0 & df$roots>1 & df$rhizo==0]="root specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="rhizo and nodule generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
#habitat
unique(habitat)


### alternative labeling
habitat<-rep(NA, 805)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="specialist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="specialist"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="generalist"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="generalist"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
unique(habitat)

############### calucating sums by hand
df_sums<-rowsum(df, habitat)
df_sums$rhizo
value <- df_sums$rhizo
value[c(2,3)]<-df_sums$nodule[c(2,3)]
value[c(6,7)]<-df_sums$roots[c(6,7)]
value
df_sums<-mutate(df_sums, value = value)

########df wrangling for venndriagram
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
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s")
svg(file="OTU_level_active_venn.svg",width = 4, height=4 )

windows(4,4)
ggvenn(
  x, 
  fill_color = mycols6[c(4,5,6)],
  stroke_size = 2, set_name_size = 4, text_size = 3, digits = 1, fill_alpha=.7
  #auto_scale = TRUE
 ) +
  ggtitle("Active Otus")
dev.off()


###########BARPLOT Abundance of habitat types ############
##make df for active ###
### don't include taxa that are only in 1 samples and less than 50 reads total
sample_data(ps)
# at least in 2 samples min reads is 10
df<-subset_samples(ps, Fraction=="BONCAT_Active")
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df<-prune_taxa(taxa_sums(df) > 0, df)
df
# 806 taxa
# rm others
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))
df
n<-row.names(df)
n[grepl("R" , n)]="rhizo"
n[grepl("E" , n)]="roots"
n[grepl("N" , n)]="nodule"

# summ by group
df<-rowsum(df, n)
df<-as.data.frame(t(df))

df<-df[-which(row.names(df)=="Others"),]

#habitat<-rep(NA, 805)
#habitat[df$nodule==0 & df$roots==0 & df$rhizo>1]="rhizo specailist"
#habitat[df$nodule==0 & df$roots>1 & df$rhizo==0]="root specialist"
#habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
#habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
#habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="rhizo and nodule generalist"
#habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
#habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
#unique(habitat)

### alternative labeling
habitat<-rep(NA, 805)
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="1 habitat"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="1 habitat"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="1 habitat"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="2 habitats"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="2 habitats"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="2 habitats"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="all habitats"
unique(habitat)

df_habitat<-cbind(df, habitat)
df_habitat<-gather(df_habitat, "Compartment", value, 1:3 )
active<-filter(df_habitat, value>0)
active$Fraction <- "Active"
active
unique(active$habitat)
######make df for total###########

df<-subset_samples(ps, Fraction=="Total_Cells")
# at least in 2 samples min reads is 10
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df<-prune_taxa(taxa_sums(df) > 0, df)
# 1293 taxa
taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))
n<-row.names.data.frame(df)
n[grepl("N" , n)]="nodule"
n[grepl("E" , n)]="roots"
n[grepl("R" , n)]="rhizo"

# summ by group
df<-rowsum(df, n)
dim(df)
df<-as.data.frame(t(df))

df<-df[-which(row.names(df)=="Others"),]

#habitat<-rep(NA, length(df$nodule))
#habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="rhizo specailist"
#habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="root specialist"
#habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="nodule specialist"
#habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="rhizo and root generalist"
#habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="rhizo and nodule generalist"
#habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="plant generalist"
#habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="hyper generalist"
#unique(habitat)

### alternative labeling
habitat<-rep(NA, length(df$nodule))
habitat[df$nodule==0 & df$roots==0 & df$rhizo>0]="1 habitat"
habitat[df$nodule==0 & df$roots>0 & df$rhizo==0]="1 habitat"
habitat[df$nodule>0 & df$roots==0 & df$rhizo==0]="1 habitat"
habitat[df$nodule==0 & df$roots>0 & df$rhizo>0]="2 habitats"
habitat[df$nodule>0 & df$roots==0 & df$rhizo>0]="2 habitats"
habitat[df$nodule>0 & df$roots>0 & df$rhizo==0]="2 habitats"
habitat[df$nodule>0 & df$roots>0 & df$rhizo>0]="all habitats"
unique(habitat)

#df wrangling for barplot ##
df_habitat<-cbind(df, habitat)
df_habitat<-gather(df_habitat, "Compartment", value, 1:3 )
total<-filter(df_habitat, value>0)
total$Fraction <- "Total"
# rm zeros???
total
# put dfs together and set factors 
total
active
df1<-rbind(active, total)
df1

#df1$habitat<-factor(df1$habitat, levels = c( "specialist",  "generalist" , "hyper generalist" ))
unique(levels(as.factor(df1$habitat)))
df1$habitat<- as.factor((df1$habitat))
df1$Compartment<-factor(df1$Compartment, levels = c( "rhizo",  "roots", "nodule" ))
unique(df1$Compartment)
df1$value <- df1$value +1 


windows(6,8)
df1%>% 
  ggplot()+
  geom_boxplot(mapping=aes(x = habitat, y = log10(value), fill=Fraction))+
  #geom_jitter(mapping=aes(x = reorder(habitat, +rhizo), y = log10(rhizo)))+
  scale_fill_manual(values= mycols6[c(1,4)])+
  facet_wrap(~Compartment, nrow = 3,
             ncol = 1)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=0.9))+
  ylab("relative abundance")+
  xlab("presence across multiple habitats")


df_habitat%>% filter(roots>1)%>%
  ggplot()+
  geom_boxplot(mapping=aes(x = reorder(habitat, +roots), y = log10(roots)), fill=mycols6[5])+
  #geom_jitter(mapping=aes(x = reorder(habitat, +roots), y = log10(roots)))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=0.9))+
  ylab("log 10 relative abundance in root")+
  xlab("habitat use group")+
  ggtitle("Roots - Active")


df_habitat%>% filter(nodule>1)%>%
  ggplot()+
  geom_boxplot(mapping=aes(x = reorder(habitat, +nodule), y = log10(nodule)), fill=mycols6[6], alpha=.7)+
  geom_jitter(mapping=aes(x = reorder(habitat, +nodule), y = log10(nodule)))+
  
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=0.9))+
  ylab("log 10 relative abundance in root")

######venn diagram of endo total vs active####
### don't include taxa that are only in 1 samples and less than 50 reads total
sample_data(ps)
# at least in 2 samples min reads is 10
df<-subset_samples(ps, Compartment=="Roots")
df<-ps_prune(df, min.samples = 2, min.reads = 10)
df<-prune_taxa(taxa_sums(df) > 0, df)
df
endo<- df # place holder
# prune to in at least 2 samples, and at least 10 reads
# 121 taxa

taxon<-tax_table(df)
df<-as.data.frame((otu_table(df)))

dim(df)
df
n<-row.names.data.frame(df)
n[grep("POS" , n)]="Active"
n[grep("SYBR" , n)]="Total"
# summ by group
df<-rowsum(df, n)

df<-t(df)
df<-as.data.frame(df)
dim(df)
head(df)

################ df wrangling for venndriagram
df[df>1] <- 1
head(df)
Active<-rownames(df[df$Active==1,])
Total<-rownames(df[df$Total==1,])
#rhizo<-rownames(df[df$rhizo==1,])
df
x <- list(
  Active = Active, 
  Total = Total
)
x

################################## venn diagram #
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
# otu level
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/venndiagram")
svg(file="OTU_level_active_venn.svg",width = 6, height=6 )

windows(6,6)
ggvenn(
  x, 
  fill_color = c("#cd622e","#4499f5"),
  stroke_size = 2, set_name_size = 5, text_size = 6, digits = 1, fill_alpha=.8
  #auto_scale = TRUE
) +
  ggtitle("Otus with at least 10 reads")
dev.off()

####who are they?
total_only<-rownames(df[df$Active==0 & df$Total==1,])
Active_only<-rownames(df[df$Active==1 & df$Total==0,])
Active_only<-as.data.frame(taxon[match(Active_only, row.names(taxon)),])
Active_only$Genus[which(Active_only$Genus==" g__"  )] <- Active_only$Family[which(Active_only$Genus==" g__" )]
Active_only$Genus[which(Active_only$Genus==""  )] <- Active_only$Family[which(Active_only$Genus=="" )]
Active_only$Genus[which(Active_only$Genus==" f__"  )] <- Active_only$Order[which(Active_only$Genus==" f__" )]
Active_only$Genus[which(Active_only$Genus==""  )] <- Active_only$Order[which(Active_only$Genus=="" )]
Active_only$Genus[which(Active_only$Genus==""  )] <- Active_only$Class[which(Active_only$Genus=="" )]
Active_only$Genus
# what is the abundance of these otus?
taxon<-tax_table(endo)
df<-as.data.frame((otu_table(endo)))

n<-row.names.data.frame(df)
n[grep("POS" , n)]="Active"
n[grep("SYBR" , n)]="Total"
# summ by group
df<-rowsum(df, n)
rowSums(df)
df<-t(df)
df<-as.data.frame(df)

df1<-as.data.frame(df[match(row.names(Active_only), row.names(taxon)),])
df<-cbind(df1, Active_only)
df<-df[order(df$Active, decreasing = TRUE),]
df%>% select(Genus, Active)

# what is the abundance of these otus only in total?
taxon<-tax_table(endo)
df<-as.data.frame((otu_table(endo)))

n<-row.names.data.frame(df)
n[grep("POS" , n)]="Active"
n[grep("SYBR" , n)]="Total"
# summ by group
df<-rowsum(df, n)

rowSums(df)
df<-t(df)
df<-as.data.frame(df)


df1<-as.data.frame(taxon[match(total_only, row.names(taxon)),])
df1
#df<-cbind(df1, _only)
#df<-df[order(df$Active, decreasing = TRUE),]
df%>% select(Genus, Active)



#####PCOA plots######
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
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Fraction, metadat$Compartment))

##------make phyloseq object with rarefied data -------#
dim(asvs.raw)
dim(asvs.r)
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
# 12855 taxa

######remove unassigned taxa ########
# taxa that are in the plant and are unassigned
# remove 
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

#Pcoa on rarefied asvs Data
# rm ctl
ps<-subset_samples(ps, Compartment !="ctl")
ps
sample_data(ps)
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
metadat2<-filter(metadat, Compartment!="ctl")
as.factor(metadat2$Compartment)
as.factor(metadat2$Fraction)
as.factor(metadat2$compartment_BCAT)
levels(as.factor(metadat2$compartment_BCAT))
unique(levels(as.factor(metadat2$Fraction)))

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s")

svg(file="Pcoa_plantvsoil_raw.svg",width = 5, height=5 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 5, height = 5)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="All Compartments",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, col=c("black"),
       pch=c(25, circle, diamond, triangle)[as.factor(metadat$Compartment)],
       lwd=1,cex=2,
       bg=c( "#cd622e", "#4499f5", "#b46db3" )[as.factor(metadat$Fraction)])

legend("topleft", legend=c( "Active"   ,   "Total_Cells", "TotalDNA"  ),
         #pch=c(circle, diamond, triangle, circle, diamond, triangle, 25, diamond),
         fill= c("#cd622e", "#4499f5", "#b46db3" ),
          title = "Fraction",
         bty = "n")

legend("top", legend=c("Nodule", "Root", "Rhizosphere", "Bulk soil"  ),
       pch=c(1,  2, 5, 6),
       #fill= c("#cd622e", "#4499f5", "#b46db3" ) ,
       title = "Compartment",       bty = "n")
dev.off()
#change pch values so they are nicer

windows(4,4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of Bray Curtis dissimilarities",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
legend("top", inset=c(-2,0), legend=c( "Active"   ,   "Total_Cells", "TotalDNA"  ),
      #pch=c(circle, diamond, triangle, circle, diamond, triangle, 25, diamond),
       fill= c("#cd622e", "#4499f5", "#b46db3" ) ,
       bty = "n")

legend("right", legend=c("Nodule", "Root", "Rhizosphere", "Bulk soil"  ),
       pch=c(circle,  triangle, diamond, 25),
       #fill= c("#cd622e", "#4499f5", "#b46db3" ) ,
       bty = "n")



ordiellipse(otus.pcoa, metadat$Compartment,  
            kind = "sd", conf=0.95, label=T, 
            #draw = "polygon",
            lwd=2, col="grey")
#
bg=c("#f0ac7d", "#cd622e", "#993203","#8fcafd" ,"#4499f5", "#0c62af","#b46db3", "#b46db3" )
dev.off()
#bulk soil = upside down tri 25
nodule = circle
root = tri
rhizo = diamond
levels(as.factor(metadat$compartment_BCAT))
as.factor(metadat$Compartment)
dev.off()

######soil####
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
metadat2<-as.data.frame(sample_data(ps2))
#metadat2<-metadat%>% filter(Compartment !=  "Nodule" & Compartment != "Roots" & Compartment!="ctl" & Fraction != "Total_DNA")

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/pcoa/soil_raw.svg",width = 4, height=4 )
windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Rhizosphere",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=(diamond),
       lwd=1,cex=2,
       bg=c(fill= c("#cd622e", "#4499f5" ) )[as.factor(metadat2$Fraction)])

ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= c("#cd622e", "#4499f5"),
            alpha = 50)
dev.off()
#ordiellipse(otus.pcoa, metadat2$Fraction,  
 #           kind = "ehull", conf=0.9, label=T, 
  #          #draw = "polygon",
   #         lwd=2, col="black")

legend("center",legend=c( "Rhizosphere BONCAT_Active", "Rhizosphere Total Cells"), 
       pch=c(circle, triangle),
       cex=1.1, 
       fill=c(blue,"#739AFF"),
       bty = "n")

legend("topright", inset=c(-0.2,0), legend=c("Rhizosphere BONCAT_Active", "Rhizosphere Total Cells"), pch=c(circle, triangle))



dev.off()

####### active verse total RDA





######nodule#####
ps
ps2<-subset_samples(ps, Compartment == "Nodule")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
df<-as.data.frame(otu_table(ps2))
sample_data(ps2)
#  146 taXA
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

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/nodule.svg",width = 4, height=4 )
windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Nodule",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=c(circle),
       lwd=1,cex=2,
       bg=c(fill= c("#cd622e", "#4499f5" ) )[as.factor(metadat2$Fraction)])

ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= c("#cd622e", "#4499f5"),
            alpha = 50)


dev.off()

####### rooots #######
ps
ps2<-subset_samples(ps, Compartment == "Roots")
ps2<-prune_taxa(taxa_sums(ps2) > 0, ps2)
any(taxa_sums(ps2) == 0)
ps2
sample_names(ps2)
df<-as.data.frame(otu_table(ps2))
#  618  taXA
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

# subset metadata
metadat2<- sample_data(ps2)

levels(as.factor(metadat2$compartment_BCAT))
levels(as.factor(metadat2$Fraction))

square <- 22
diamond <- 23
triangle <- 24
circle <- 21

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/pcoa/roots.svg",width = 4, height=4 )
windows(title="PCoA on asvs- Bray Curtis", width = 4, height = 4)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="Roots",xlab=paste("PCoA1(",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p[,1:2],
       pch=triangle,
       lwd=1,cex=2,
       bg=c(fill= c("#cd622e", "#4499f5" ) )[as.factor(metadat2$Fraction)])

ordiellipse(otus.pcoa, metadat2$Fraction,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            #lwd=.1,
            col= c("#cd622e", "#4499f5"),
            alpha = 50)


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
# Permanova just plant #
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




#####HEATMAP START#########

#####--------heatmap  family level activity-------########

##select active taxa 
# at least in 2 samples min reads is 50
df<-subset_samples(ps, Fraction=="BONCAT_Active")
df<-ps_prune(df, min.samples = 2, min.reads = 100)
df # 510 taxa
taxon<-as.data.frame(tax_table(df))
df<-as.data.frame(t(as.data.frame(otu_table(df))))
head(df)
# 717 taxa
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

#rm other columns
taxainfo <- select(df, Phyla, Domain, Class,  Order, Family )
taxainfo<-unique(taxainfo)
df<-select(df, -Phyla, -Domain, -Class, -Order, -Genus, -Species, -otu)

#aggregate
df<-aggregate(cbind(C10E.POS_S60,  C10N.POS_S65, C1E.POS_S31, C1N.POS_S61,  C2E.POS_S32,  C2N.POS_S62,  C5E.POS_S33,  C5N.POS_S63, C10R.POS_S30, C1R.POS_S27, C2R.POS_S28, C5R.POS_S29 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 113 families

#### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()
df$root <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)

# for now rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>%   glimpse()
# call the family column famliy
df<-mutate(df, Family=otu) %>% select(., -otu)
head(df)
dim(df)
#### I'm little curious of the the taxonomy of these families
taxainfo<-left_join(df, taxainfo)
dim(taxainfo)
head(taxainfo)

### plan:
#### for taxa that have at least 1000 reads
### identified to the order level but not the famliy level will be added to the tree manually


#importtree
#tree <- read.tree("2022.10.phylogeny.asv.nwk")
#tree <- read.tree("2022.10.phylogeny.id.nwk")
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/trees/GTDB")

tree = read.tree("family.nwk")
tree
# 1057 tips
length(tree$tip.label) # look at the tip labels 
# modify tip labels
tree$tip.label <- paste0(" f__", tree$tip.label)

length(intersect(unique(df$Family), tree$tip.label)) # Apply setdiff function to see what's missing from the tree
# 58 intersect
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

# make matrix

m <- df
summary(m)
m<-log10(m)
m[m== "-Inf"] <- 0
m[m==0] <- NA
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m

# heat map
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/figures/16s/heatmap")
#svg(file="Family_level_active_heatmap.svg",width = 8, height=8 )

# brenden's function
#col_dendro = as.dendrogram(hclust(dist(t(m))))
#windows(10,10)
#        heatmap.phylo(x = m, Rowp = tree.short, Colp = as.phylo(as.hclust(col_dendro)))

# ggtree
windows(6,6)
p<- ggtree(tree.short, branch.length = .001) + 
 geom_tiplab(size=3, align=TRUE, linesize=0, offset = -.2) + 
 theme_tree2()+
 xlim_tree(1) 
gheatmap(p, m, 
         colnames=FALSE,
         legend_title="active taxa", offset = .5) +
  scale_x_ggtree() + 
  scale_fill_gradient(low= "#fffaa2", high =  "#bb0000", aesthetics = "fill", na.value = "white",
                      name="Abundance in Active")+
  ggtitle("Families in Active Community")
dev.off()

### no phylogeny ###########
#make matrix
m <- df
m<-log10(m)
m[m== "-Inf"] <- 0
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m
lighhyellow == "#fffaa2"
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "#fcdb47", "#bb0000"))(n = 199)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.01,length=2),  # for red
               seq(0.1,0.8,length=98),           # for yellow
               seq(0.81,5.1,length=100))             # for green


# clustering agrorithm
distance = dist(m, method = "euclidean")
distance
cluster = hclust(distance, method = "ward.D2")
cluster$labels
# useful https://sebastianraschka.com/Articles/heatmaps_in_r.html
dim(m)
row.names(m)

windows(6,12)
svg(filename = "heatmap_active.svg", width = 6, height = 12)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "row",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          labRow = cluster$labels,
          margins = c(7,7),
          cexRow= .65,
          cexCol = 1,
          #key.title = "log 10 Abundance",
          #key.xlab = "")
          key = FALSE) 
  #title("Active", line= -4)
dev.off()  
  
 

#####--------heatmap family level total-------########
##select active taxa 
# at least in 2 samples min reads is 50
sample_data(ps)
ps1<-subset_samples(ps, Fraction=="Total_Cells")
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 100)
any(taxa_sums(ps1) == 0)
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))
dim(df)
# 596 taxa

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
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/trees/GTDB")

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
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s/")

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


### no phylogeny ###########
#make matrix
m <- df
m<-log10(m)
m[m== "-Inf"] <- 0
m<-m[,c(3,2,1)]
m<-as.matrix(m)
m

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("white", "#8fcafd", "#0c328a"))(n = 199)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0.01,length=2),  # for red
               seq(0.1,0.8,length=98),           # for yellow
               seq(0.81,5.1,length=100))             # for green


# clustering agrorithm
distance = dist(m, method = "euclidean")
distance
cluster = hclust(distance, method = "ward.D2")
cluster$labels
# useful https://sebastianraschka.com/Articles/heatmaps_in_r.html


windows(6,12)
svg(filename = "heatmap_total.svg", width = 6, height = 12)
heatmap.2(m, 
          col = my_palette,
          breaks = col_breaks,
          density.info="none",
          dendrogram = "row",
          trace="none",         # turns off trace lines inside the heat ma
          Rowv = as.dendrogram(cluster),# turns off density plot inside color legend
          Colv="NA",
          labRow = cluster$labels,
          margins = c(7,7),
          cexRow= .65,
          cexCol = 1,
          #key.title = "log 10 Abundance",
          #key.xlab = "")
          key = FALSE) 
#title("Active", line= -4)
dev.off()  










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

##Import tree##
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

##### otu heatmap active ####  
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


############################## heat map habitat use type active ######

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
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>% filter(otu!="Other") %>%glimpse()
df<-mutate(df, Family=otu) %>% select(., -otu)
summary(df)
head(df)


habitat<-rep(NA, 138)
habitat[df$nodule==0 & df$root==0 & df$rhizo>0]="rhizospecialist"
habitat[df$nodule==0 & df$root>0 & df$rhizo==0]="rootspecialist"
habitat[df$nodule>0 & df$root==0 & df$rhizo==0]="nodulespecialist"
habitat[df$nodule==0 & df$root>0 & df$rhizo>0]="rhizoandrootgeneralist"
habitat[df$nodule>0 & df$root>0 & df$rhizo==0]="plantgeneralist"
habitat[df$nodule>0 & df$root>0 & df$rhizo>0]="hypergeneralist"
habitat[df$nodule>0 & df$root==0 & df$rhizo>0]="rhizoandnodulegeneralist"

df$habitat <- habitat
df_habitat <-df
head(df_habitat)

#who are they?
n<-df%>% filter(habitat=="rootspecialist") %>% select(Family) 
n
n<-df%>% filter(habitat=="rhizospecialist") %>% select(Family) 
n
n<-df%>% filter(habitat=="rhizoandrootgeneralist") %>% select(Family) 
n
n<-df%>% filter(habitat=="plantgeneralist") %>% select(Family) 
n
n<-df%>% filter(habitat=="hypergeneralist") %>% select(Family) 
n

#make a little tree with the habitat types
df_habitat<-df_habitat[order(df_habitat$rhizo, df$root),]
df<-df_habitat%>% select(habitat,Family)

colnames(df) <- c("x", "y")
df<-df[order(df$x),]
row.names(df)<- NULL

newick<-df2newick(df, innerlabel = TRUE)
tree <- read.tree(text=newick)
tree

###get the correct order of the tree
###remove the space
df_habitat$Family<-sub("\\ ", "", df_habitat$Family)
target<-tree$tip.label
df_habitat<-df_habitat[match(target, df_habitat$Family),]
row.names(df_habitat) <- df_habitat$Family

############### make matrix
m<-df_habitat
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


#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/heatmap2.svg", width = 10, height=10 )
# And make the plot with phylogeny
#windows(8,8)
#heatmap.phylo(x = m, Rowp = tree, Colp = as.phylo(as.hclust(col_dendro)))+
#  dev.off()

tree.c<-read.tree(text =  "((rhizo:2):1, (root:2):1, (nodule:2):1);")
head(m)
dim(m)
tree.c
colnames(m)
tree$tip.label
row.names(m)

heatmap(m, Rowv=tree.short, Colv=tree.c)

heatmap.phylo(x = m, Rowp = tree, Colp = tr)
dev.off
# no phylogeny
#make a dendrogram              
row_dendro = as.dendrogram(hclust(dist((m))))
windows(10,10)
heatmap.phylo(x = m, Rowp = as.phylo(as.hclust(row_dendro)), Colp = as.phylo(as.hclust(col_dendro)))
#legend(legend(x="right", legend=c("min", "med", "max"),fill=
dev.off

#### this is cool we make this with just the top ten taxa in each place
n1<-df_habitat %>% filter(habitat=="rhizospecialist") %>%
  .[order(.$rhizo, decreasing = TRUE),] %>% .[1:10,] %>% .$Family

n2<-df_habitat %>%filter(habitat=="hypergeneralist")%>% .$Family

n3 <-  df_habitat %>% filter(habitat=="rhizoandrootgeneralist") %>%
  .[order(.$rhizo, decreasing = TRUE),] %>% .[1:10,] %>% .$Family

n4 <-  df_habitat %>% filter(habitat=="rootspecialist") %>%
  .[order(.$root, decreasing = TRUE),] %>% .[1:10,] %>% .$Family

keep<-c(n1, n2, n3, n4)

#### little heatmap by habitat use typ

#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, keep) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE)
#
tree.short
target<-tree.short$tip.label
df_habitat<-df_habitat[match(target, df_habitat$Family),]
row.names(df_habitat) <- df_habitat$Family
dim(df_habitat)

##make matrix
m<-df_habitat
m<-as.matrix(m[,c(3,2,1)])
dim(m)
#matrix needs a log transform
m<-m+1
m<-log10(m)
#m[m== "-Inf"] <- 0
#m[m== 0] <-NA
head(m)
#summary(m)
#heatmap(m)

#plot
windows(14,8)
heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)

##### heatmap by habitat use total #######

##select active taxa 
# at least in 2 samples min reads is 50
ps1<-subset_samples(ps, Fraction=="Total_Cells")
ps1<-ps_prune(ps1, min.samples = 2, min.reads = 50)
df<-as.data.frame(t(as.data.frame(otu_table(ps1))))
taxon<-as.data.frame(tax_table(ps))
dim(df)
# 971 taxa
####### aggregate to family level again
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

#aggregate
df<-aggregate(cbind(C10N.SYBR_S26, C10R.SYBR_S20, C1E.SYBR_S21,  C1N.SYBR_S13,  C1R.SYBR_S16,  C2E.SYBR_S22,
                    C2N.SYBR_S15,  C2R.SYBR_S17 ,  C5E.SYBR_S23,  C5R.SYBR_S18,  C7E.SYBR_S24,  C7N.SYBR_S25,  C7R.SYBR_S19 ) ~ Family, data = df, FUN = sum, na.rm = TRUE)
row.names(df) <- df$Family
dim(df) # 145 families

### summarize by compartment 
df$rhizo <-   rowMeans(df %>% dplyr::select(contains("R.SYB")))
df$nodule <-   rowMeans(df %>% dplyr::select(contains("N.SYB")))
df$root <-   rowMeans(df %>% dplyr::select(contains("E.SYB")))
df<-df%>% select(c(nodule, root, rhizo)) %>% mutate(otu= row.names(df))
head(df)
# for now rm the unknown family row
df<-df%>%filter(otu!='') %>% filter(otu!=" f__")  %>% filter(otu!="Other") %>%glimpse()
df<-mutate(df, Family=otu) %>% select(., -otu)
summary(df)
head(df)
dim(df)


habitat<-rep(NA, 143)
habitat[df$nodule==0 & df$root==0 & df$rhizo>0]="rhizospecialist"
habitat[df$nodule==0 & df$root>0 & df$rhizo==0]="rootspecialist"
habitat[df$nodule>0 & df$root==0 & df$rhizo==0]="nodulespecialist"
habitat[df$nodule==0 & df$root>0 & df$rhizo>0]="rhizoandrootgeneralist"
habitat[df$nodule>0 & df$root>0 & df$rhizo==0]="plantgeneralist"
habitat[df$nodule>0 & df$root>0 & df$rhizo>0]="hypergeneralist"
habitat[df$nodule>0 & df$root==0 & df$rhizo>0]="rhizoandnodulegeneralist"

df$habitat <- habitat
df_habitat <-df
head(df_habitat)

#who are they?
n<-df%>% filter(habitat=="rootspecialist") %>% select(Family) 
n
n<-df%>% filter(habitat=="rhizospecialist") %>% select(Family) 
n
n<-df%>% filter(habitat=="rhizoandrootgeneralist") %>% select(Family) 
n
n<-df%>% filter(habitat=="plantgeneralist") %>% select(Family) 
n
n<-df%>% filter(habitat=="hypergeneralist") %>% select(Family) 
n

# make tree
#make a little tree with the habitat types
df_habitat<-df_habitat[order(df_habitat$rhizo, df$root),]
df<-df_habitat%>% select(habitat,Family)

colnames(df) <- c("x", "y")
df<-df[order(df$x),]
row.names(df)<- NULL

newick<-df2newick(df, innerlabel = TRUE)
tree <- read.tree(text=newick)
tree

###get the correct order of the tree
###remove the space
df_habitat$Family<-sub("\\ ", "", df_habitat$Family)
target<-tree$tip.label
df_habitat<-df_habitat[match(target, df_habitat$Family),]
row.names(df_habitat) <- df_habitat$Family

############### make matrix
m<-df_habitat
m<-as.matrix(m[,c(3,2,1)])
dim(m)
#log transform
m<-m+1
m<-log10(m)
#m[m== "-Inf"] <- 0
#m[m== 0] <-NA
head(m)
summary(m)
heatmap(m)

#make dendrogram
tree.c<-read.tree(text =  "((rhizo:2):1, (root:2):1, (nodule:2):1);")


heatmap.phylo(x = m, Rowp = tree, Colp = tr)


# make plot #
# load in the function for making a heatmap with the tree #
heatmap.phylo <- function(x, Rowp, Colp, ...) {
  l = length(seq(.1, 6, .1))
  pal = colorRampPalette(c("#8fcafd", "#0c328a"))(l)
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

# plot
heatmap.phylo(x = m, Rowp = tree, Colp = tr)


# shorten
#### this is cool we make this with just the top ten taxa in each place
n1<-df_habitat %>% filter(habitat=="rhizospecialist") %>%
  .[order(.$rhizo, decreasing = TRUE),] %>% .[1:10,] %>% .$Family

n2<-df_habitat %>%filter(habitat=="hypergeneralist")%>% .$Family

n3 <-  df_habitat %>% filter(habitat=="rhizoandrootgeneralist") %>%
  .[order(.$rhizo, decreasing = TRUE),] %>% .[1:10,] %>% .$Family

n4 <-  df_habitat %>% filter(habitat=="rootspecialist") %>% .$Family

n5<-df_habitat %>%filter(habitat=="rhizoandnodulegeneralist")%>% .$Family
unique(df_habitat$habitat)

keep<-c(n1, n2, n3, n4, n5)
keep
#### little heatmap by habitat use typ

#shorten phylogeny to match what is in our data frame
asvs_remove<-setdiff(tree$tip.label, keep) #asvs we don't want
tree.short<-drop.tip(tree, asvs_remove) # remove asvs we don't need
plot(tree.short, no.margin=TRUE)
#
tree.short
target<-tree.short$tip.label
df_habitat<-df_habitat[match(target, df_habitat$Family),]
row.names(df_habitat) <- df_habitat$Family
dim(df_habitat)

##make matrix
m<-df_habitat
m<-as.matrix(m[,c(3,2,1)])
dim(m)
#matrix needs a log transform
m<-m+1
m<-log10(m)
#m[m== "-Inf"] <- 0
#m[m== 0] <-NA
head(m)
#summary(m)
#heatmap(m)

#plot
windows(14,8)
heatmap.phylo(x = m, Rowp = tree.short, Colp = tr)

### I'm not loving this plot make use gg tree instead. 
windows(8,8)
p <- ggtree(tree.short, branch.length = 'none') + 
  #xlim_tree(2) +
  geom_tiplab(size=3, align=FALSE, linesize=1, offset = 0) + 
  theme_tree2()
p
gheatmap(p, m, 
         colnames=FALSE, legend_title="total taxa") +
  scale_x_ggtree() + 
  # scale_color_gradient(l
  scale_fill_gradient(low="#8fcafd" , high = "#0c328a", aesthetics = "fill", na.value = "white",
                      name="Abundance in TOTAL")+
  ggtitle("Habitat use")
dev.off()








################### phylogeny graveyard #################

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


#changing my data's names so they match the tree #

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




####  endosphyte taxa exploration ########

# remove taxa that do not have a phyla and are not in the soil
soil <- subset_samples( ps, Compartment=="Rhizosphere" | Compartment=="Bulk_soil" | Compartment="ctl")
soil
# which taxa are zero in the soil?
notinsoil<-prune_taxa(taxa_sums(soil) ==0, soil)
notinsoil
# 5168 taxa not in soil
# okay for not in soil, who is is not bacteria?
#<_subset_taxa(notinsoil, Domain!="d__Bacteria" )
# 34 taxa
####### remove taxa with no phyla

# okay for not in soil, who esn't have a ohyla
remove<-subset_taxa(ps, Phyla=="")
remove
# 482 taxa
#tax_table(remove)

########## remove these guys

#remove<-otu_table(remove)
badtaxa<-taxa_names(remove)

alltaxa<-taxa_names(ps)
mytaxa <- alltaxa[!(alltaxa %in% badtaxa)]

ps<-prune_taxa(mytaxa, ps )
ps
# 12373 taxa if we just remove the ones that are not in soil
# 11838 taxa

###### remake figure


# grab data
taxon<- as.data.frame(tax_table(ps))
df<-t(otu_table(ps))
df<-cbind(df, taxon)

# summarize by phyla
df<-aggregate(cbind(BEADS_S67 ,    C10B.DNA_S4  , C10E.POS_S60  ,C10N.POS_S65 , C10N.SYBR_S26 , C10R.DNA_S12 ,  C10R.POS_S30 ,
                    C10R.SYBR_S20 , C1E.POS_S31,   C1E.SYBR_S21,  C1N.POS_S61,   C1N.SYBR_S13,  C1R.DNA_S8,    C1R.POS_S27,  
                    C1R.SYBR_S16,  C2B.DNA_S1,    C2E.POS_S32,   C2E.SYBR_S22,  C2N.POS_S62,   C2N.SYBR_S15,  C2R.DNA_S10,  
                    C2R.POS_S28,   C2R.SYBR_S17,  C5B.DNA_S2,    C5E.POS_S33,   C5E.SYBR_S23,  C5N.POS_S63,   C5R.DNA_S11,  
                    C5R.POS_S29 ,   C5R.SYBR_S18,  C7B.DNA_S3,    C7E.POS_S34,   C7E.SYBR_S24,  C7N.POS_S64,   C7N.SYBR_S25, 
                    C7R.DNA_S14,   C7R.SYBR_S19,  CTL_S66,       S10.DNA_S9,    S2.DNA_S5,     S3.DNA_S6,     S8.DNA_S7
) ~ Phyla, data = df, FUN = sum, na.rm = TRUE)

head(df)
df$Phyla[c(1,2)] <- "other"

# again so there is only 1 other
df<-aggregate(cbind(BEADS_S67 ,    C10B.DNA_S4  , C10E.POS_S60  ,C10N.POS_S65 , C10N.SYBR_S26 , C10R.DNA_S12 ,  C10R.POS_S30 ,
                    C10R.SYBR_S20 , C1E.POS_S31,   C1E.SYBR_S21,  C1N.POS_S61,   C1N.SYBR_S13,  C1R.DNA_S8,    C1R.POS_S27,  
                    C1R.SYBR_S16,  C2B.DNA_S1,    C2E.POS_S32,   C2E.SYBR_S22,  C2N.POS_S62,   C2N.SYBR_S15,  C2R.DNA_S10,  
                    C2R.POS_S28,   C2R.SYBR_S17,  C5B.DNA_S2,    C5E.POS_S33,   C5E.SYBR_S23,  C5N.POS_S63,   C5R.DNA_S11,  
                    C5R.POS_S29 ,   C5R.SYBR_S18,  C7B.DNA_S3,    C7E.POS_S34,   C7E.SYBR_S24,  C7N.POS_S64,   C7N.SYBR_S25, 
                    C7R.DNA_S14,   C7R.SYBR_S19,  CTL_S66,       S10.DNA_S9,    S2.DNA_S5,     S3.DNA_S6,     S8.DNA_S7
) ~ Phyla, data = df, FUN = sum, na.rm = TRUE)


# ?? look at how we plot the other one?
df<-gather(df, "sample", value, 2:43 )

#remove zeros
# anything that is less than 1% equals other
df<-df[df$value!=0,]
df
#df1<-df # place holder
#df1$Phyla[df1$value<1] <- "other"

## who are the other?
#<-df[grep("E.", df$sample),]

df1<-aggregate(cbind(value) ~ sample+Phyla, data = df1, FUN = sum, na.rm =TRUE)

windows(10,10)
df%>% 
  ggplot(aes(fill=Phyla, y=value, x=sample)) + 
  geom_bar(position="fill", stat= "identity")+
  #scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  scale_fill_viridis(discrete = TRUE) +
  ggtitle("Top phyla") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

dev.off()




###########binomial model #########

###make data frame##
#presence in plant ~ abundance in rhizosphere
# at least in 2 samples min reads is 10
#df<-ps_prune(df, min.samples = 2, min.reads = 10)
#df<-prune_taxa(taxa_sums(df) > 0, df)
#dim(df)
# 1293 taxa
# don't include taxa that are super rare less then 5 read
#ps
#sample_data(ps)
ps1
ps1<-subset_samples(ps, Compartment != "ctl"& Fraction != "Total_DNA")
#ps1 12802 taxa
taxon<- as.data.frame(tax_table(ps1))
df<-as.data.frame(otu_table(ps1))

#### seperate out ecach rep
#               abundance_active  abundance_total present_in_plant
#otu #1 rep 1
#otu #1 rep 2
n<-row.names.data.frame(df)
n[grepl("C10" , n)]="5"
n[grepl("C1" , n)]="1"
n[grepl("C2" , n)]="2"
n[grepl("C5" , n)]="3"
n[grepl("C7" , n)]="4"
df$rep <- n

#rep 1
df1<-filter(df, rep=="1")
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
inplant[df1$nodule>0 | df1$roots>0 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep1")
colnames(df1) <- c("Active", "Total", "inplant")
rep1<-df1
rep1

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
inplant[df1$nodule>0 | df1$roots>0 ]<-"1"
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
inplant[df1$nodule>0 | df1$roots>0 ]<-"1"
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
inplant[df1$nodule>0 | df1$roots>0 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep4")
colnames(df1) <- c("Total", "inplant")
df1
rep4<-df1

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
inplant[df1$nodule>0 | df1$roots>0 ]<-"1"
inplant[df1$nodule==0 & df1$roots==0 ]<-"0"
df1$inplant <- inplant
df1<-df1 %>% select(c(-nodule, -roots))
row.names(df1)<-paste0(row.names(df1), "Rep5")
#change column names
colnames(df1) <- c("Active", "Total", "inplant")
rep5<-df1

### put them together
df1<-rbind(rep1, rep2) %>% rbind(., rep3) %>% rbind(., rep5)
df1$inplant <- as.numeric(df1$inplant)
colnames(df1)
dim(df1)
# 51208 3
######SCATTERPLOT by inplant is red ######
head(df1)

#plot
filter(df1, Total >10 & Active >10) %>% 
  ggplot(aes(x=log10(Total) , y=log10(Active), col = as.factor(inplant))) + 
  geom_jitter(stat="identity", position="identity") +
  scale_color_manual(values = c("black", "red")) +
  geom_abline(slope=1)+
  theme_bw()+
  labs(x="Log 10 Abundance in Total Viable Cell population",y="log 10 Abundance in Active population")


# remove taxa that have less than 10 reads
df2<-df1[df1$Active>10 & df1>10,]
dim(df2)
#4706 10 both
#16348 rows
#42651 present in both

#logistic regression
#total
#plot 1:
ggplot(df2, aes(y=inplant, x=Total))+
 geom_point()
#Model 1: 
fit <- glm(as.numeric(df2$inplant) ~ df2$Total, family = binomial)
summary(fit)
# AIC 1262.6
ggplot()+
  geom_smooth(data = df2, aes(x = Active, y = inplant),
              method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  geom_point(data = df2, aes(x = Active, y = inplant))

# both

#Model 2: 
fit <- glm(as.numeric(df2$inplant) ~ df2$Total + df2$Active + df2$Active*df2$Total, family = binomial)
summary(fit)
# active is a way stronger predictor
# AIC 1231.4
ggplot()+
  geom_smooth(data = df2, aes(x = Total, y = inplant),
              method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  geom_point(data = df2, aes(x = Total, y = inplant))

#active
#plot 3:
ggplot(df2, aes(y=inplant, x=Active))+
 geom_point()
#Model 3: 
fit <- glm(as.numeric(df2$inplant) ~ df2$Active, family = binomial)
summary(fit)
# AIC 1235.7
ggplot()+
  geom_smooth(data = df2, aes(x = Active, y = inplant),
              method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  geom_point(data = df2, aes(x = Active, y = inplant))




######LFC from ancom######
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/16s")

tab<-read_delim("table_ancom.txt")

tab<-tab%>% select(!contains("DNA"))

##otus table
df<-as.data.frame(t(otu_table(ps)))
head(df)
df$otu<-row.names(df)
# taxon table
taxon<-as.data.frame(tax_table(ps1))
df<-left_join(taxon, df)
df<-left_join(tab, df)
dim(df)
# negative lower in total than active 
# that confusing 
# multiplyu by -1 to make negative is lower in total
df<-mutate(df, LFC_Fraction_Total_cells_Active<-LFC_Fraction_Total_cells_Active*-1)
df

## pool columns with plant
colnames(df)
i<-rowSums(df%>%select(contains("E.") | contains("N.")))
i[i>0] <-1
i
df$inplant <- i

df<-df%>%select(otu, inplant, LFC_Fraction_Total_cells_Active,  pval_Fraction_Total_cells_Active, Phyla, Family, Genus)
df<-mutate(df, LFC=LFC_Fraction_Total_cells_Active, pval=pval_Fraction_Total_cells_Active) %>% select(-LFC_Fraction_Total_cells_Active)


#Plot 1:
p<-ggplot(df, aes(y=inplant, x=LFC))
p + geom_point()


#Plot 2:
df%>% filter() %>%
ggplot(aes(y=pval, x=LFC, col=as.factor(inplant)))+
  scale_color_manual(values = c("grey", "red"))+
  geom_point()+
  theme_bw(base_size = 16)


# Plot 3:
ggplot()+
  geom_smooth(data = df, aes(x = LFC, y = inplant),
              method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  geom_point(data = df, aes(x = LFC, y = inplant))+
  theme_bw(base_size = 16)
  
  
#Plot 4 
  df%>% filter(LFC>0) %>%
  ggplot(aes(y=pval, x=LFC, col=as.factor(inplant)))+
  scale_color_manual(values = c("grey", "red"))+
  geom_point()+
  theme_bw(base_size = 16)+

# plot 5

ggplot()+
  geom_smooth(data = df%>% filter(LFC>0), aes(x = LFC, y = inplant),
              method = "glm", method.args = list(family = "binomial"), se = FALSE)+
  geom_point(data = df%>% filter(LFC>0), aes(x = LFC, y = inplant))+
  theme_bw(base_size = 16)

  
  #Model 4:
  fit <- glm(as.numeric(df$inplant) ~ df$LFC, family = binomial)
  summary(fit)
