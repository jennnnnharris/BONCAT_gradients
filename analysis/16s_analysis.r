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
## Convert OTU numbers to percentages ##
#otus.perc<-otus.r/rowSums(otus.r)*100

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
# 14593 taxa

otus<-as.data.frame(t(as.data.frame(otu_table(ps))))
taxon<-as.data.frame(tax_table(ps))

#################------ calculating inactive fraction---------------####
#### remove taxa from bulk soil
ps1<- subset_samples(ps,Fraction !="Total_DNA" )
ps1<-prune_taxa(taxa_sums(ps1) > 0, ps1)
any(taxa_sums(ps1) == 0)
ps1

# 7185 taxa

# subset total cell and active fraction
ps.Total<- subset_samples(ps1,Fraction =="Total_Cells" )
otu_total<-as.data.frame(t(as.data.frame(otu_table(ps.Total))))
# boncat pos table doesn't have C7Rpos, so I'll remove that from the table

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

#sum<-rowSums(otu_active)
# length 7185
# taxa1 10000
# taxa2 0
# taxa3 0 


#subset
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
#
#
colnames(otu_inactive)
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

#just checking the row names and col names look right
rownames(otus)
colnames(otus)
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
  ggplot(aes(x=reorder(compartment_BCAT,-Observed), y=Observed, fill = Fraction))+
  geom_boxplot() +
  scale_fill_manual(values = c("grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab("Number of ASVs")+
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
  ggplot(aes(x=compartment_BCAT, y=Observed, fill = Fraction, col= Fraction))+
  geom_boxplot() +
  scale_colour_manual(values = c( "orange",  "black", "black"))+
  scale_fill_manual( values = c("gold", "grey27", "lightgrey"))+
  geom_jitter(width = .1, size=1 )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  ylab("Number of ASVs")+
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
svg(file="figures/16s/pcoa/pcoa_all_raw.svg",width = 7, height=6 )

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




# pc 3 and 4
# Store coordinates for axes 3 adn 4 in new variable #
otus.p <- otus.pcoa$points[,1:4]

# subset
plant.p<-otus.p[which(metadat$Fraction != "Total_Cells"),]
plant.p<-plant.p[which(metadat_in$Compartment !="Rhizosphere" & metadat_in$Compartment !="Bulk_Soil" ),]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]
pe3<-perc.exp[3]
pe4<-perc.exp[4]
pe5<-perc.exp[5]
#to make coloring things easier I'm gong to added a combined fractionXboncat column no sure if i need this
metadat_pl<-metadat%>% filter(Compartment !="Rhizosphere") %>% filter(Compartment!="Bulk_Soil") %>% filter(Fraction != "Total_Cells")

as.factor(metadat_pl$Fraction)
as.factor(metadat_pl$Compartment)
as.factor(metadat_pl$compartment_BCAT)
## Plot ordination with factors coloured and shaped as you like # 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/Pc1.pc2.plants.svg",width = 6, height=6 )
windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(1,2), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA1 (",round(pe1, 2),"% variance explained)"),
         ylab=paste("PCoA2 (",round(pe2,2),"% variance explained)"))
points(plant.p,
       pch=c(22,21,0,24)[as.factor(metadat_pl$Fraction)],
       lwd=1,cex=2,
       bg=c("black", "white", "#DC267F", "#DC267F",  "#FFB000", "#FFB000")[as.factor(metadat_pl$compartment_BCAT)])
# roots = gold
# nodules = magenta
# ctl = white
#active = circle
# inactive = trangle

legend("top",legend=c("Flow cyto control",  " PCR control",  "NoduleBONCAT_Active" , "NoduleInactive",  "RootsBONCAT_Active",  "RootsInactive"), 
       pch=c(15,0,1,2,1,2 ),
       cex=1.1, 
       col=c("black", "black", "#DC267F", "#DC267F",  "#FFB000", "#FFB000"  ),
       inset = c(.2, 0),
       bty = "n")

dev.off()

#############--------- plants pc3 pc 4--------------------

## Plot ordination with factors coloured and shaped as you like # 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
#svg(file="figures/16s/Pc3.pc4.plants.svg",width = 7, height=6 )

windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(3,4), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA3 (",round(pe3, 2),"% variance explained)"),
         ylab=paste("PCoA4 (",round(pe4,2),"% variance explained)"))
points(plant.p[,3:4],
       pch=c(22,21,0,24)[as.factor(metadat_pl$Fraction)],
       lwd=1,cex=2,
       bg=c("black", "white", "#DC267F", "#DC267F",  "#FFB000", "#FFB000")[as.factor(metadat_pl$compartment_BCAT)])
# roots = gold
# nodules = magenta
# ctl = white
#active = circle
# inactive = trangle

legend("top",legend=c("Flow cyto control",  " PCR control",  "NoduleBONCAT_Active" , "NoduleInactive",  "RootsBONCAT_Active",  "RootsInactive"), 
       pch=c(15,0,1,2,1,2 ),
       cex=1.1, 
       col=c("black", "black", "#DC267F", "#DC267F",  "#FFB000", "#FFB000"  ),
       inset = c(.2, 0),
       bty = "n")

dev.off()

##################### soil and rhizosphere #############
# Store coordinates for axes 3 adn 4 in new variable #
otus.p <- otus.pcoa$points[,1:4]

# subset
soil.p<-otus.p[which(metadat$Fraction != "Total_Cells"),]
soil.p<-soil.p[which(metadat_in$Compartment !="Nodule" & metadat_in$Compartment !="Roots" ),]

# Calculate % variance explained by each axis #
otus.eig<-otus.pcoa$eig
perc.exp<-otus.eig/(sum(otus.eig))*100
pe1<-perc.exp[1]
pe2<-perc.exp[2]
pe3<-perc.exp[3]
pe4<-perc.exp[4]
pe5<-perc.exp[5]
#to make coloring things easier I'm gong to added a combined fractionXboncat column no sure if i need this
metadat_so<-metadat%>% filter(Compartment !="Nodule") %>% filter(Compartment!="Roots") %>% filter(Fraction != "Total_Cells")

as.factor(metadat_so$Fraction)
as.factor(metadat_so$Compartment)
as.factor(metadat_so$compartment_BCAT)
## Plot ordination with factors coloured and shaped as you like # 
#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/Pc2.pc3.soil.svg",width = 6, height=6 )
windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(2,3), type="none", main="PCoA of Bray Curtis",xlab=paste("PCoA2(",round(pe2, 2),"% variance explained)"),
         ylab=paste("PCoA3 (",round(pe3,2),"% variance explained)"))
points(soil.p[,2:3],
       pch=c(22,21,0,24, 23)[as.factor(metadat_so$Fraction)],
       lwd=1,cex=2,
       bg=c("black","#739AFF", "white", "#785EF0", "#785EF0",  "#785EF0")[as.factor(metadat_so$compartment_BCAT)])
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

#setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/")
svg(file="figures/16s/Pc3.pc4.soil.svg",width = 6, height=6 )
windows(title="PCoA on plant asvs- Bray Curtis", width = 7, height = 6)
ordiplot(otus.pcoa,choices=c(3,4), type="none", main="PCoA of 16S OTU Bray Curtis",xlab=paste("PCoA3(",round(pe3, 2),"% variance explained)"),
         ylab=paste("PCoA4 (",round(pe4,2),"% variance explained)"))
points(soil.p[,3:4],
       pch=c(22,21,0,24, 23)[as.factor(metadat_so$Fraction)],
       lwd=1,cex=2,
       bg=c("black","#739AFF", "white", "#785EF0", "#785EF0",  "#785EF0")[as.factor(metadat_so$compartment_BCAT)])
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
