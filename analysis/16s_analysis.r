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


####### OTU 97% data##################
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

####------data wrangleing and importing into Phyloseq-----------#####
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

## Determine minimum available reads per sample ##
min(rowSums(otu.t))
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
ps <- phyloseq(Workshop_taxo, Workshop_OTU, Workshop_metadat)

##  take a look at PS object
print(ps)
# we a get a Plyloseq object with  7727 taxa

#######------ remove chloroplasts-------######
# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that arent in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
# 7516 taxa
# output df 
otus<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))

########------diversity figure---------########
# set wd for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

# diversity 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

svg(file="figures/16s/diversity.svg",width = 10, height=4 )
#windows()
plot_richness(ps, "Fraction", measures = c("Observed","Shannon", "Simpson", "InvSimpson")) +
  geom_boxplot(aes(fill = "Fraction")) + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578"))
dev.off()

# diversity of active microbes in each fraction
rich<-cbind(rich, metadat)
rich<-as.data.frame(rich)
colnames(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))

reorder(compartment_BCAT, -Observed)

svg(file="figures/16s/total_alpha_diversity.svg",width = 5, height=4 )
windows()
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", Fraction!="BONCAT_Active")%>%
  ggplot(aes(x=reorder(compartment_BCAT,-Observed), y=Shannon, fill = Fraction))+
  geom_boxplot() +
  scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylab()+
  xlab("Compartment")
dev.off()

rich<- rich %>%  filter(Plant!="NOPLANT", Fraction!="Inactive")

rich$compartment_BCAT <-factor(rich$compartment_BCAT, levels = c("Bulk_SoilTotal_DNA", "RhizosphereTotal_DNA", "RhizosphereTotal_Cells", "RhizosphereBONCAT_Active",
                                                                 "RootsTotal_Cells"   ,  "RootsBONCAT_Active" , "NoduleTotal_Cells" ,  "NoduleBONCAT_Active" ))          

# total + active 
svg(file="figures/16s/active_alpha_diversity_2.svg",width = 6, height=4 )
windows()
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=Shannon, fill = Fraction, col= Fraction))+
  geom_boxplot() +
  scale_colour_manual(values = c( "orange",  "black", "black"))+
  scale_fill_manual( values = c("gold", "grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylab("Number of ASVs")+
  xlab("Compartment")
dev.off()


# look at what's in the endophyte

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

##------- log fold change from active to inactive ----------
## calculated a total fraction 
## remove taxa from bulk soil
ps1<- subset_samples(ps,Fraction !="Total_DNA"& Fraction!="beads" & Fraction !="ctl" ) 
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
# 3128 taxa

#Who are the most abundant taxa minus the total DNA
# mayeb should evaluted top taxa in active ??? idk???
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
tt<-as.data.table(otu_table(Top_tax))
hist(tt$f18c54e493d5d5ad534b9c3100216416)

# subset total cell and active fraction
ps.Total<- subset_samples(ps1,Fraction =="Total_Cells" )
otu_total<-as.data.frame(t(as.data.frame(otu_table(ps.Total))))
# boncat pos table doesn't have C7Rpos, so I'll remove that from the table
head(otu_total)
#remove this sample because I don't have a boncat sample for it
otu_total<-select(otu_total, -C7R.SYBR_S19)
#dim(otu_total)
# length #3128
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
#dim(otu_active)
#sum<-rowSums(otu_active)
# length 3128
# taxa1 10000
# taxa2 0
# taxa3 0 

# log fold change from active to inactive 
# add 1 to everything (absent in active & absent in total = no change)

otu_log2<-log2(otu_active+1/(otu_total+1))

n<-c("C10N", "C10R" ,"C1E" , "C1N" , "C1R" , "C2E" , "C2N" , "C2R" , "C5E" , "C5R" , "C7E" , "C7N" )
colnames(otu_log2)<-n
head(otu_log2)
# check distribution
otu_log2
hist(otu_log2$C10N) # most things didn't change because most taxa not present 
hist(otu_log2$C10R)
hist(otu_log2$C1E, breaks = 20)
# insert otu columns
otus<-row.names(otu_log2)
otu_log2<-mutate(otu_log2, otus=otus)

#####wrangling log2 data #####
#transpose
#t_log2_otu<- t(otu_log2)
#import metadata for summarising 
#metadat_l<- read.delim("16s/metadata_log2.txt", sep="\t", header = T, check.names=FALSE)
#y<-colnames(otu_log2)
#rownames(metadat_l) <- y
#metadat_l<-as.data.frame(metadat_l)
#insert metadata to data frame
#t_log2_otu<-cbind(metadat_l, t_log2)
#sumarise by rep
#t_log2_otu_short<- t_log2_otu %>% group_by(Fraction) %>% summarise_if(is.numeric,median) %>% select(-REP)

# filter dataset for top taxa
dim(tt)
#I'm curious who these asvs are so I'm going to add there names in
tt200<-as.data.frame(t(tt))
asvs<-row.names(tt200)
tt200<-mutate(tt200, asvs=asvs)
#filter the log2 dataset by the top taxa
otu_log2_200<-inner_join(otu_log2, tt200, by= "asvs")
tt200$asvs
# inner join keeps to many columns
otu_log2_200<-otu_log2_200[,(1:12)]
dim(otu_log2_200)
otu_log2_200<-as.matrix(otu_log2_200)
row.names(otu_log2_200)<-tt200$asvs

#### heatmap of top taxa ############
# load in the function for making a heatmap with the tree #
heatmap.phylo <- function(x, Rowp, Colp, ...) {
  l = length(seq(-4.9, 5, 0.1))
  pal = colorRampPalette(c('#2166ac', '#f7f7f7', '#b2182b'))(l)
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
        xaxs="i", yaxs="i", axes=FALSE, xlab="",ylab="",breaks=seq(-5,5,0.1))
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", yaxs="i", xlim=c(0,2), ylim=yl)
  text(rep(0,nrow(x)),1:nrow(x), row_order, pos=4, family='Helvetica',
       cex=1, xpd=NA)
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", xaxs="i", ylim=c(0,2), xlim=xl)
  text(1:ncol(x),rep(2,ncol(x)), col_order, srt=90, adj=c(1,.5), family='Helvetica',
       cex=1.5)
}


## Read in the tree file made by jenn
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")
tree = read.tree("otu_output/tree.nwk")

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
tt<-as.data.table(otu_table(Top_tax))
hist(tt$f18c54e493d5d5ad534b9c3100216416)

# filter dataset for top taxa
dim(tt)

#filter the log2 dataset by the top taxa
otu_log2_200<-inner_join(otu_log2, tt200, by= "asvs")
tt200$asvs
# inner join keeps to many columns
otu_log2_200<-otu_log2_200[,(1:12)]
dim(otu_log2_200)
otu_log2_200<-as.matrix(otu_log2_200)
row.names(otu_log2_200)<-tt200$asvs

# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(otu_log2_200)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
# tree is so big. filter and make it smaller
asvs200 # asvs we want
asvs<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(asvs, asvs200) #asvs we don't want
tree.all<-drop.tip(tree, asvs_remove) # remove asvs we don't need

# And )make the plot....
#pdf(file="../figures/Fig2_heatmap24wk_Phylo.pdf",width = 5,height=8, useDingbats=FALSE)
windows(6,10)
heatmap.phylo(x = m, Rowp = tree, Colp = as.phylo(as.hclust(col_dendro)))
#dev.off()

##### well this heat maps isn't super helpful
##### because the Nodules and roots are so much less diverse than the soil

# i am curious what taxa is enriched in the rhizosphere I should just focus on rhizosphere
# if i'm curious about what makes a good endo
# make is this figure for the taxa that are enrich in the plant
# and for the taxa that are enriched in the soil
# if taxa that are really enriched in the plant are also really enriched in the soil 
# then it could show that soil fitness helps colonize the plant?
# but if it's super different (enriched in the plant almost not present in the soil)
# then that could suport that there are plant specialists (who are they?)

##top 10 in nodule
#subset
#What;s in the nodule?
nod<-subset_samples(ps, Compartment=="Nodule")
nod<-prune_taxa(taxa_sums(nod) > 0, nod)
any(taxa_sums(nod) == 0)
nod
# only 91 taxa
topN=10
most_abundant_taxa = sort(taxa_sums(nod), TRUE)[1:topN]
Top_tax = prune_taxa(names(most_abundant_taxa),nod)

length(get_taxa_unique(Top_tax, "Class"))
length(get_taxa_unique(Top_tax, "Genus"))
nod_tt<-as.data.table(otu_table(Top_tax))
#subset for these 100 taxa
dim(nod_tt)
nod_tt
# transpose df
tt100<-as.data.frame(t(nod_tt))
# make otu column to join by in short df
otus100<-row.names(tt100)
tt100<-mutate(tt100, otus=otus100)
#join df
otu_log2_100<-inner_join(otu_log2, tt100, by= "otus")
otu_log2_100<-otu_log2_100[,1:12]
n<-colnames(otu_log2_100)
# make a matrix
otu_log2_100<-as.matrix(otu_log2_100)
row.names(otu_log2_100)<-tt100$otus
row.names(otu_log2_100)
# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(otu_log2_100)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
# tree is so big. filter and make it smaller
otus100<- tt100$otus
# asvs we want
otus<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(otus, otus100) #asvs we don't want
tree.nod<-drop.tip(tree, asvs_remove) # remove asvs we don't need

# And make the plot....
#pdf(file="../figures/Fig2_heatmap24wk_Phylo.pdf",width = 5,height=8, useDingbats=FALSE)
windows(6,10)
heatmap.phylo(x = m, Rowp = tree.nod, Colp = as.phylo(as.hclust(col_dendro)))
#dev.off()



##top 100 in enodphyte
#subset
#What;s in the endoule?
endo<-subset_samples(ps, Compartment=="Roots")
endo<-prune_taxa(taxa_sums(endo) > 0, endo)
any(taxa_sums(endo) == 0)
endo
# only 406 taxa
topN=100
most_abundant_taxa = sort(taxa_sums(endo), TRUE)[1:topN]
Top_tax = prune_taxa(names(most_abundant_taxa),endo)

length(get_taxa_unique(Top_tax, "Class"))
length(get_taxa_unique(Top_tax, "Genus"))
endo_tt<-as.data.table(otu_table(Top_tax))
#subset for these 100 taxa
dim(endo_tt)
endo_tt
# transpose df
tt100<-as.data.frame(t(endo_tt))
# make otu column to join by in short df
otus100<-row.names(tt100)
tt100<-mutate(tt100, otus=otus100)
#join df
otu_log2_100<-inner_join(otu_log2, tt100, by= "otus")
otu_log2_100<-otu_log2_100[,1:12]
n<-colnames(otu_log2_100)
# make a matrix
otu_log2_100<-as.matrix(otu_log2_100)
row.names(otu_log2_100)<-tt100$otus
row.names(otu_log2_100)
otu_log2_100
# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(otu_log2_100)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
# tree is so big. filter and make it smaller
otus100<- tt100$otus
# asvs we want
otus<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(otus, otus100) #asvs we don't want
length(otus)
length(otus100)
tree.endo<-drop.tip(tree, asvs_remove) # remove asvs we don't need

# And make the plot....
#pdf(file="../figures/Fig2_heatmap24wk_Phylo.pdf",width = 5,height=8, useDingbats=FALSE)
windows(6,10)
heatmap.phylo(x = m, Rowp = tree.endo, Colp = as.phylo(as.hclust(col_dendro)))
#dev.off()





















##top 100 in soil
#subset
#What's in the soil?
rhiz<-subset_samples(ps, Compartment=="Rhizosphere")
rhiz<-prune_taxa(taxa_sums(rhiz) > 0, rhiz)
any(taxa_sums(rhiz) == 0)
rhiz
# 4682 taxa
topN=100
most_abundant_taxa = sort(taxa_sums(rhiz), TRUE)[1:topN]
Top_tax = prune_taxa(names(most_abundant_taxa),rhiz)
length(get_taxa_unique(Top_tax, "Class"))
length(get_taxa_unique(Top_tax, "Genus"))
rhiz_tt<-as.data.table(otu_table(Top_tax))
#subset for these 100 taxa
dim(rhiz_tt)
rhiz_tt
# transpose df
tt100<-as.data.frame(t(rhiz_tt))

# make otu column to join by in short df
otus100<-row.names(tt100)
tt100<-mutate(tt100, otus=otus100)
#join df
otu_log2_100<-inner_join(otu_log2, tt100, by= "otus")
otu_log2_100<-otu_log2_100[,1:12]
n<-colnames(otu_log2_100)
# make a matrix
otu_log2_100<-as.matrix(otu_log2_100)
row.names(otu_log2_100)<-tt100$otus
row.names(otu_log2_100)
otu_log2_100
# Create the matrix and get the column dendrogram for the heatmap from it.
m<-structure(otu_log2_100)
#make a dendrogram              
col_dendro = as.dendrogram(hclust(dist(t(m))))
# tree is so big. filter and make it smaller
otus100<- tt100$otus
# asvs we want
otus<-row.names(otu.raw) # all the asvs
asvs_remove<-setdiff(otus, otus100) #asvs we don't want
length(otus)
length(otus100)
tree.rhiz<-drop.tip(tree, asvs_remove) # remove asvs we don't need

# And make the plot....
#pdf(file="../figures/Fig2_heatmap24wk_Phylo.pdf",width = 5,height=8, useDingbats=FALSE)
windows(6,10)
heatmap.phylo(x = m, Rowp = tree.rhiz, Colp = as.phylo(as.hclust(col_dendro)))
#dev.off()
# um who are these taxa tho?


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

















######-------------Asvs-----------------##############
## Set the working directory; modify to your own ###
setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")

### Import Data ###
taxon <- read.table("16s/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("16s/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("16s/metadata.txt", sep="\t", header = T, check.names=FALSE)

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
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

# take a look at PS object
print(ps)
# we a get a Plyloseq object with  15027 taxa

## Determine minimum available reads per sample ##
min(rowSums(otus.t))

### Rarefy to obtain even numbers of reads by sample ###
set.seed(336)
otus.r<-rrarefy(otus.t, 41610)
############------------ rarefaction curve------------############

#rarecurve
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


######--- recode metadata----- ########
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))

#####------make phyloseq object with rarefied data -------#####

otus.phyloseq<- (otus.r)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otus.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_OTU <- otu_table(otus.phyloseq, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

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

otus<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))


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

###@# heatmap of top taxa

#library (gplots)
library(ape)
## Read in the tree file made by jenn
tree = read.tree("16s/tree.nwk")
# my tree just has ASV names think


# load in the function for making a heatmap with the tree #
heatmap.phylo <- function(x, Rowp, Colp, ...) {
  l = length(seq(-4.9, 5, 0.1))
  pal = colorRampPalette(c('#2166ac', '#f7f7f7', '#b2182b'))(l)
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
        xaxs="i", yaxs="i", axes=FALSE, xlab="",ylab="",breaks=seq(-5,5,0.1))
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", yaxs="i", xlim=c(0,2), ylim=yl)
  text(rep(0,nrow(x)),1:nrow(x), row_order, pos=4, family='Helvetica',
       cex=1, xpd=NA)
  par(mar=rep(0,4))
  plot(NA, axes=FALSE, ylab="", xlab="", xaxs="i", ylim=c(0,2), xlim=xl)
  text(1:ncol(x),rep(2,ncol(x)), col_order, srt=90, adj=c(1,.5), family='Helvetica',
       cex=1.5)
}

# Create the matrix and get the column dendrogram for the heatmap from it.
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












##---- calculating inactive fraction-------
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


#####------make phyloseq object with rarefied data + inactive fraction -------#####

otus.phyloseq<- t(otus)
#head(otus.phyloseq)
metadat<-as.matrix(metadat)
y<-colnames(otus)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_OTU <- otu_table(as.matrix(otus.phyloseq), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon)) # this taxon file is from the prev phyloseq object length = 14833
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)


#test it worked
sample_names(ps)
print(ps)
# 14593 taxa

#####1. Calculate diversity#######

# set wd for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

# diversity 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

svg(file="figures/16s/diversity.svg",width = 10, height=4 )
#windows()
plot_richness(ps, "Fraction", measures = c("Observed","Shannon", "Simpson", "InvSimpson")) +
geom_boxplot(aes(fill = "Fraction")) + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578"))
dev.off()

# diversity of active microbes in each fraction
rich<-cbind(rich, metadat)
rich<-as.data.frame(rich)
colnames(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))

reorder(compartment_BCAT, -Observed)

svg(file="figures/16s/total_alpha_diversity.svg",width = 5, height=4 )
windows()
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive", Fraction!="BONCAT_Active")%>%
  ggplot(aes(x=reorder(compartment_BCAT,-Observed), y=Shannon, fill = Fraction))+
  geom_boxplot() +
  scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylab()+
  xlab("Compartment")
dev.off()

rich<- rich %>%  filter(Plant!="NOPLANT", Fraction!="Inactive")

rich$compartment_BCAT <-factor(rich$compartment_BCAT, levels = c("Bulk_SoilTotal_DNA", "RhizosphereTotal_DNA", "RhizosphereTotal_Cells", "RhizosphereBONCAT_Active",
           "RootsTotal_Cells"   ,  "RootsBONCAT_Active" , "NoduleTotal_Cells" ,  "NoduleBONCAT_Active" ))          

# total + active 
svg(file="figures/16s/active_alpha_diversity_2.svg",width = 6, height=4 )
windows()
rich %>%
  filter(Plant!="NOPLANT", Fraction!="Inactive")%>%
  ggplot(aes(x=compartment_BCAT, y=InvSimpson, fill = Fraction, col= Fraction))+
  geom_boxplot() +
  scale_colour_manual(values = c( "orange",  "black", "black"))+
  scale_fill_manual( values = c("gold", "grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylab("Number of ASVs")+
  xlab("Compartment")
dev.off()


# look at what's in the endophyte

rich %>%
  filter(BONCAT!="DNA", Fraction!="ctl", Compartment== "Roots", REP!=1)%>%
  ggplot(aes(x=Fraction, y=Observed, col=REP))+
  #geom_boxplot() +
  #geom_line( x=Fraction, y=Observed, col= rep )
  #scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  ylab("Number of ASVs")+
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

dim(filter(total, sum>1))

dim(filter(root_active, sum>1))



#phylogenic diversity
#https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#alpha-diversity-calculations


#####2. Phyloseq and manipulation by taxonomy ####

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

#------- 3. Run PCoA analysis of entire dataset ------


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

