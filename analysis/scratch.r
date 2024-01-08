
# scratch 


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
#df1<-df # place holder
# order columns

# make really low abundance taxa other
df1$Phyla[df1$value<1] <- "other"
df1<-aggregate(cbind(value) ~ sample+Phyla, data = df1, FUN = sum, na.rm =TRUE)
head(df1)
df1<-df1[order(df1$sample),]
head(df1)

#how many phyla?
length(unique(df1$Phyla))
# 18 families 

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
# n cells
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


