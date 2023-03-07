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
otus.raw <- read.table("16s/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("16s/metadata.txt", sep="\t", header = T, check.names=FALSE)

## Transpose OTU table ##
otus.t <- t(otus.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order otu table
otus.t<-otus.t[order(row.names(otus.t)),]

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

## Convert OTU numbers to percentages ##
#otus.perc<-otus.r/rowSums(otus.r)*100

#set colors
#mycols=c("grey27","grey","#9795ff","#4406e2", "#ffb4f6","#e20a8f", "#695e00", "#6eda0a", "#0a7416" )

# rhizo, endo, nod, 
#green , grey, purple , grey, pink, grey
#mycols3= c("#3aaf04", "grey","#4406e2", "grey" , "#ff50c8", "grey")

# recode some names so they are easier to understand
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
#metadat$Compartment<-factor(metadat$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))
#to make coloring things easier I'm gong to added a combined fractionXboncat column 
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))
metadat$compartment_BCAT<-as.factor(metadat$compartment_BCAT)

#####------make phyloseq object with rarefied data -------#####

otus.phyloseq<- t(otus.r)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otus.raw)
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
# 15027 taxa

# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")
# get rid of taxa that aren; in any samples
ps<-prune_taxa(taxa_sums(ps) > 0, ps)
any(taxa_sums(ps) == 0)
ps
# 14833 taxa

otus<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))

#################------ calculating inactive fraction---------------####
#### remove taxa from bulk soil
ps1<- subset_samples(ps,Fraction !="Total_DNA" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1
# 7289 taxa

# subset total cell and active fraction
ps.Total<- subset_samples(ps1,Fraction =="Total_Cells" )
ps.Total
otu_total<-as.data.frame(t(as.data.frame(otu_table(ps.Total))))
otu_total
# length #7289
#       C10N.SYBR_S26 C10R.SYBR_S20 C1E.SYBR_S21 C1N.SYBR_S13 
# taxa1  139581          1430        62828       126909  
# taxa2  19537           179         8553        17541
# taxa3  78                0         4336        22797

ps.Active<- subset_samples(ps1,Fraction =="BONCAT_Active" )
ps.Active
otu_active<-as.data.frame(t(as.data.frame(otu_table(ps.Active))))
otu_active<-rowSums(otu_active)
# length 7289
# taxa1 10000
# taxa2 0
# taxa3 0 

#subset
#Inactive<- total- active
# all the taxa that are in the total fraction but not in active are the "inactive fraction"
# this gives all the row that have zero values, length = 4001
inactive<-which(otu_active==0)
otu_inactive<-otu_total[inactive,]
#length 4001
#       C10N.SYBR_S26 C10R.SYBR_S20 C1E.SYBR_S21 C1N.SYBR_S13 
# taxa1  0              47          0            0  
# taxa2  0              0           0            0
# taxa3  0              22          0            0
#
#

# change colnames
n<-c("C10N.inactive", "C10R.inactive", "C1E.inactive"  ,"C1N.inactive" , "C1R.inactive",  "C2E.inactive" , "C2N.inactive"  ,"C2R.inactive" , "C5E.inactive", 
"C5R.inactive",  "C7E.inactive",  "C7N.inactive",  "C7R.inactive" )
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

rownames(otus)
colnames(otus)
# append metadata


metadat <- read.delim("16s/metadata_w_inactive.txt", sep="\t", header = T, check.names=FALSE)
metadat<-metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule', 'ctl'='ctl'))
metadat<-metadat[, c(1,3:6)]
metadat<-metadat%>% mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl'))

#to make coloring things easier I'm gong to added a combined fractionXboncat column no sure if i need this
metadat<-mutate(metadat, compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction))
metadat<-select(metadat, -BONCAT)

## order metadata
metadat<-metadat[order(metadat$SampleID),]
metadat$SampleID
## order otu table
otus<-otus[, order(colnames(otus))]
colnames(otus)

## Convert OTU numbers to percentages ##
otus.perc<-otus/rowSums(otus)*100


#####------make phyloseq object with rarefied data + inactive fraction -------#####

otus.phyloseq<- t(otus)
head(otus.phyloseq)
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

#####1. Calculate diversity#######

# set wd for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")


# diversity 
rich<-estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson", "InvSimpson" ))

svg(file="figures/16s/diversity.svg",width = 10, height=6 )
#windows()
plot_richness(ps, "Fraction", measures = c("Observed","Shannon", "Simpson", "InvSimpson")) +
geom_boxplot(aes(fill = "Fraction")) + scale_fill_manual(values = c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578"))
dev.off()

# diversity of active microbes in each fraction

rich<-cbind(rich, metadat)
rich<-as.data.frame(rich)
colnames(rich)
rich$Compartment<-factor(rich$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))


svg(file="figures/16s/alpha_diversity.svg",width = 5, height=4 )
windows()
rich %>%
  filter(BONCAT!="DNA", Fraction!="ctl")%>%
  ggplot(aes(x=compartment_BCAT, y=Observed))+
  geom_boxplot() +
  #scale_fill_manual(values = c("grey27", "lightgrey"))+
  #geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab("Number of ASVs")
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
tax_table(root_total)
total<-as.data.table(otu_table(root_total))
t(total)

root_active<-subset_samples(ps, compartment_BCAT=="RootsBONCAT_Active")
root_active<-prune_taxa(taxa_sums(root_active) > 0, root_active)
any(taxa_sums(root) == 0)
tt<-as.data.table(tax_table(root_active))
ot<-as.data.table(otu_table(root_active))
t(ot)
root_active<-cbind(ot,tt)


#phylogenic diversity
#https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#alpha-diversity-calculations


#####2. Phyloseq and manipulation by taxonomy ####

#check n taxa
rank_names(ps)

#interacting with phyloseq object
sample_variables(ps)
length(sample_variables(ps))

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
ggplot(d1, aes(OTUabundance)) + 
  geom_histogram() +
  ggtitle("Histogram of Total Counts") + 
  xlim(0, 20000) + ylim (0,50) + theme_bw()
# most taxa are rare, but some are really abundant

#select most abundant taxa
topN = 10
most_abundant_taxa = sort(taxa_sums(ps), TRUE)[1:topN]
print(most_abundant_taxa)
GP20 = prune_taxa(names(most_abundant_taxa), ps)
length(get_taxa_unique(GP20, "Class"))
print(get_taxa_unique(GP20, "Phyla"))
print(get_taxa_unique(GP20, "Family"))
tax_table(GP20)

#negative ctl

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

#beads
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

#nodule
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

#genus in nodule  
# c( " Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium",  " Ensifer" , " Pseudomonas"  " Staphylococcus") 

####100% plots ######
rhizobiaceae <- subset_taxa(ps, Family ==  " Rhizobiaceae"  )
rhizobia.sum<-rowSums(otu_table(rhizobiaceae))
psuedo <- subset_taxa(ps, Family == " Pseudomonadaceae" )
psuedo.sum<-rowSums(otu_table(psuedo))
staph <- subset_taxa(ps, Family == " Staphylococcaceae")
staph.sum<-rowSums(otu_table(staph))
burk <- subset_taxa(ps, Family == " Burkholderiaceae"  )
burk.sum<-rowSums(otu_table(burk))
sphing <- subset_taxa(ps, Family== " Sphingomonadaceae" )
sphing.sum<-rowSums(otu_table(sphing))
micro <- subset_taxa(ps, Family ==  " Micrococcaceae"   )
micro.sum<-rowSums(otu_table(micro))

other<-rowSums(otus.t)- (rhizobia.sum+psuedo.sum+staph.sum+burk.sum+sphing.sum+micro.sum)

  
phyl.mat<-cbind(rhizobia.sum,psuedo.sum,staph.sum,burk.sum,sphing.sum, micro.sum, other)
print(phyl.mat)
phyl.ag<-aggregate(phyl.mat~metadat$Compartment,FUN=mean)

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


###### import percent abundance into phyloseq #####

otus.phyloseq<- t(otus.perc)
taxon<-taxon[,1:7]
metadat<-as.matrix(metadat)
y<-colnames(otus)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_OTU <- otu_table(as.matrix(otus.perc), taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_OTU,Workshop_metadat)

#test it worked
sample_names(ps)
print(ps)

# remove chloroplast DNA
ps<-subset_taxa(ps, Class!=" Chloroplast")
ps<-subset_taxa(ps, Genus!=" Mitochondria")
ps<-subset_taxa(ps, Genus!=" Chloroplast")

#### 100% plots

rhizobiaceae <- subset_taxa(ps, Family ==  " Rhizobiaceae"  )
rhizobia.sum<-rowSums(otu_table(rhizobiaceae))
psuedo <- subset_taxa(ps, Family == " Pseudomonadaceae" )
psuedo.sum<-rowSums(otu_table(psuedo))
staph <- subset_taxa(ps, Family == " Staphylococcaceae")
staph.sum<-rowSums(otu_table(staph))
burk <- subset_taxa(ps, Family == " Burkholderiaceae"  )
burk.sum<-rowSums(otu_table(burk))
sphing <- subset_taxa(ps, Family== " Sphingomonadaceae" )
sphing.sum<-rowSums(otu_table(sphing))
micro <- subset_taxa(ps, Family ==  " Micrococcaceae"   )
micro.sum<-rowSums(otu_table(micro))

other<-100-(rhizobia.sum+psuedo.sum+staph.sum+burk.sum+sphing.sum+micro.sum)


phyl.mat<-cbind(rhizobia.sum, psuedo.sum, staph.sum, burk.sum, sphing.sum, micro.sum, other)
print(phyl.mat)
phyl.mat<-aggregate(phyl.mat~metadat$Compartment,FUN=mean)

colnames(phyl.mat)[1]<- "Compartment"
# data wrangling  
phy.df.t<-gather(phyl.mat, "taxa", value, 2:8)

phy.df.t %>% group_by(Compartment, taxa) %>%
summarise(mean = mean(value), n = n())

# Stacked + percent

svg(file="figures/16s/top_taxa.svg",width = 6, height=5 )

ggplot(phy.df.t, aes(fill=taxa, y=value, x=Compartment)) + 
  geom_bar(position="fill", stat= "identity")+
  scale_fill_manual(values= c("#26808f","#6499b5", "#9db1d3","#d1cbe9", "#d2a0d0","#dc6e9c","#d43d51")) +
  ggtitle("Top Families") +
  theme_bw(base_size = 14) +
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
otus.pcoa <- cmdscale(otus.bray, k=(ncol(otus)-1), eig=TRUE)

# Store coordinates for first two axes in new variable #
otus.p <- otus.pcoa$points[,1:2]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]

## Plot ordination with factors coloured and shaped as you like # 
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")

#svg(file="figures/16s/Pcoa.svg",width = 6, height=6 )

windows(title="PCoA on OTUs - Bray Curtis")
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(otus.p, col=c("black"),
       pch=c(20,21,22,23,24)[as.factor(metadat$Fraction)],
       lwd=1,cex=2,
       #bg=c("#003f5c","grey", "#bc5090", "#ffa600")[as.factor(metadat$BONCAT)])
       bg=c( "#52311f","grey", "#f538de","#ffb4f6", "#276602", "#5ef507", "#dbfcc7", "#4406e2", "#a483f7"  )[metadat$compartment_BCAT])

#dark brown, grey,  nod(light pink, darkpink ), rhizo(darkest green, dark green, light green,), endo (dark purple, light purple),
# BulkDNA ctlctl EndoPOS EndoSYBR NodPOS NodSYBR RhizoDNA RhizoPOS RhizoSYBR

legend("top",legend=c("Bulk_SoilTotal_DNA" ,  "ctl"  ,"NoduleBONCAT_Active" ,"NoduleTotal_Cells", 
                      "RhizosphereBONCAT_Active","RhizosphereTotal_Cells", "RhizosphereTotal_DNA" ,"RootsBONCAT_Active"   ,   "RootsTotal_Cells"),
       pch=c(16,16,15,15,18, 18,17,17,17),
       cex=1.1, 
       col=c( "#52311f","grey", "#f538de","#ffb4f6", "#276602", "#5ef507", "#dbfcc7", "#4406e2", "#a483f7"  ))

#dev.off()

#------------active verse inactive---------

# subset
p<-otus.p[which(metadat$Fraction != "Total_Cells"),]

# subset metadata
metadat_in<-metadat%>% filter(Fraction !="Total_Cells") 
as.factor(metadat_in$Compartment)
as.factor(metadat_in$Fraction)
as.factor(metadat_in$compartment_BCAT)
levels(as.factor(metadat_in$compartment_BCAT))

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/Pcoa.svg",width = 6, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis")
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA1 (",round(pe1,2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(p, col=c("black"),
       pch=c(21,22,23,24, 25)[as.factor(metadat_in$Fraction)],
       lwd=1,cex=2,
       bg=c("#739AFF", "white", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000")[as.factor(metadat_in$compartment_BCAT)])


legend("top",legend=c( "Bulk_SoilTotal_DNA" ,      "ctlctl"   ,  "NoduleBONCAT_Active"    ,  "NoduleInactive"     ,      "RhizosphereBONCAT_Active",
                       "RhizosphereInactive","RhizosphereTotal_DNA" ,    "RootsBONCAT_Active"   ,    "RootsInactive"),
       pch=c(2,0,1, 5, 1, 5 , 2, 1, 5  ),
       col= c("#739AFF", "black", "#DC267F", "#DC267F", "#785EF0", "#785EF0" , "#785EF0", "#FFB000", "#FFB000")       )

dev.off()



#------------subset by groups------------
# Store coordinates for axes 3 adn 4 in new variable #
otus.p <- otus.pcoa$points[,1:4]

# subset
plant.p<-otus.p[which(metadat$Fraction != "Total_Cells"),]
plant.p<-plant.p[which(metadat_in$Compartment == "Nodule"| metadat_in$Compartment == "Roots" ),]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe2<-perc.exp[2]
pe3<-perc.exp[3]
pe4<-perc.exp[4]
pe5<-perc.exp[5]
#to make coloring things easier I'm gong to added a combined fractionXboncat column no sure if i need this
metadat_pl<-metadat%>% filter(Compartment !="Rhizosphere") %>% filter(Compartment!="Bulk_Soil") %>% filter(Fraction != "Total_Cells")

as.factor(metadat_pl$Fraction)
as.factor(metadat_pl$Compartment)
## Plot ordination with factors coloured and shaped as you like # 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/Pcoa.svg",width = 6, height=6 )
windows(title="PCoA on plant asvs- Bray Curtis")
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA1 (",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(plant.p, display = Sample_ID
       pch=c(21,22,22)[as.factor(metadat_pl$Fraction)],
       lwd=1,cex=2,
       col=c("black", "#DC267F", "#FFB000")[as.factor(metadat_pl$Compartment)]
       text(display = metadat_pl$SampleID))

# roots = gold
# nodules = magenta
# ctl = white
#active = circle
# inactive = square

legend("top",legend=c("Bulk_SoilTotal_DNA" ,  "ctl"  ,"NoduleBONCAT_Active" ,"NoduleTotal_Cells", 
                      "RhizosphereBONCAT_Active","RhizosphereTotal_Cells", "RhizosphereTotal_DNA" ,"RootsBONCAT_Active"   ,   "RootsTotal_Cells"),
       pch=c(16,16,15,15,18, 18,17,17,17),
       cex=1.1, 
       col=c( "#52311f","grey", "#f538de","#ffb4f6", "#276602", "#5ef507", "#dbfcc7", "#4406e2", "#a483f7"  ))

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
  
