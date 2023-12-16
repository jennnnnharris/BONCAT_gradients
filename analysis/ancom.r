#Jenn
# sep 2023

#library(BiocManager)
#BiocManager::install("microbiome")

#knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
#                      fig.width = 6.25, fig.height = 5)

#install.packages("remotes")
#remotes::install_github("Russel88/MicEco")
#install.packages("MicEco")

library(ANCOMBC)
library(MicEco)

library(tidyverse)
library(microbiome)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

##########example data #######

data(atlas1006, package = "microbiome")
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(atlas1006)
tse
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

tse

out = ancombc(data = tse, assay_name = "counts", 
              tax_level = "Family", phyloseq = NULL, 
              formula = "age + region + bmi", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global
# ancombc also supports importing data in phyloseq format
# tse_alt = agglomerateByRank(tse, "Family")
# pseq = makePhyloseqFromTreeSummarizedExperiment(tse_alt)
# out = ancombc(data = NULL, assay_name = NULL,
#               tax_level = "Family", phyloseq = pseq,
#               formula = "age + region + bmi",
#               p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
#               group = "region", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
#               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
#               n_cl = 1, verbose = TRUE)

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "Age", "NE - CE", "SE - CE", 
             "US - CE", "Overweight - Obese", "Lean - Obese")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
# standard error
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
#test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
# adjusted p value
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")

# data viz
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, age, ageSE) %>%
  dplyr::filter(age != 0) %>% 
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(direct = ifelse(age > 0, "Positive LFC", "Negative LFC"))
df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id)
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

p_age = ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = age, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = age - ageSE, ymax = age + ageSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_age
# data viz obese verse lean
df_fig_bmi = df_lfc %>% 
  filter(bmioverweight != 0 | bmilean != 0) %>%
  transmute(taxon_id, 
            `Overweight vs. Obese` = round(bmioverweight, 2),
            `Lean vs. Obese` = round(bmilean, 2)) %>%
  pivot_longer(cols = `Overweight vs. Obese`:`Lean vs. Obese`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2
p_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_bmi


###################import my data #############


## Set the working directory; ###
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/16s")


### Import Data ###
taxon <- read.table("asv_level_output/greengenes/taxonomy.txt", sep="\t", header=T, row.names=1)
asvs.raw <- read.table("asv_level_output/greengenes/feature-table.tsv", sep="\t", header=T, row.names = 1 )
metadat <- read.delim("metadata.txt", sep="\t", header = T, check.names=FALSE)

## Transpose ASVS table ##
asvs <- t(asvs.raw)
## order metadata
metadat<-metadat[order(metadat$SampleID),]
## order asvs table
asvs<-asvs[order(row.names(asvs)),]

###--- recode metadata----- #
metadat<- metadat%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule')) %>%
  mutate(Fraction=recode(BONCAT, 'DNA'= 'Total_DNA', 'SYBR'= 'Total_Cells', 'POS'='BONCAT_Active', 'ctl'= 'ctl')) %>%
  select(-BONCAT) %>%
  select(-SampleID) %>%
  mutate(compartment_BCAT = paste0(metadat$Compartment, metadat$Fraction)) %>%
  glimpse()



#prep for phyloseq
asvs.phyloseq<- (asvs)
taxon<-taxon[,1:7]
taxon$otu<-rownames(taxon)
metadat<-as.matrix(metadat)
y<-colnames(asvs.raw)
rownames(metadat) <- y
metadat<-as.data.frame(metadat)

#import it phyloseq
Workshop_ASVS <- otu_table(asvs, taxa_are_rows = FALSE)
Workshop_metadat <- sample_data(metadat)
Workshop_taxo <- tax_table(as.matrix(taxon))
ps <- phyloseq(Workshop_taxo, Workshop_ASVS,Workshop_metadat)

#test it worked
#sample_names(ps)
print(ps)
# 12855 taxa


#####################data frame prep sum and mean ########
ps1



df<-as.data.frame(t(otu_table(ps1)))
# BAR CHART with TOP 50 OTUS
#nodule
df$nodule.total.sum <-   rowSums(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
df$nodule.xbcat.sum <-   rowSums(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()

# roots
df$root.total.sum <-   rowSums(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
df$root.xbcat.sum <-   rowSums(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()

#rhizo
df$rhizo.total.sum <-   rowSums(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df$rhizo.bcat.sum <-   rowSums(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$rhizo.DNA.sum <-   rowSums(df %>% dplyr::select(contains("R.DNA"))) %>%   glimpse()
colnames(df)

#put means too
#nodule
df$nodule.total.mean <-   rowMeans(df %>% dplyr::select(contains("N.SYBR"))) %>%   glimpse()
df$nodule.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("N.POS"))) %>%   glimpse()

# roots
df$root.total.mean <-   rowMeans(df %>% dplyr::select(contains("E.SYBR"))) %>% glimpse()
df$root.xbcat.mean <-   rowMeans(df %>% dplyr::select(contains("E.POS"))) %>%   glimpse()

#rhizo
df$rhizo.total.mean <-   rowMeans(df %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df$rhizo.bcat.mean <-   rowMeans(df %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df$rhizo.DNA.mean <-   rowMeans(df %>% dplyr::select(contains("R.DNA"))) %>%   glimpse()
colnames(df)




########ANCOM #########
sample_data(ps)
ps1<-subset_samples(ps, Compartment=="Rhizosphere")
ps1<-ps_prune(ps1, min.samples = 3, min.reads = 10)
sample_data(ps1)

any(taxa_sums(ps1) == 0)

ps1
#1823 taxa if you remove things that are only present in at least 3 samples and just rhizo

#taxa if you remove things that are only present in 1 sample
#5995 if you want to rm taxa that are super rare
sample_data(ps1)

ps2 = mia::makeTreeSummarizedExperimentFromPhyloseq(ps1)



out1 = ancombc(data = ps2, assay_name = "counts", 
              tax_level = "otu", phyloseq = NULL, 
              formula = "Fraction", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "Fraction", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)

res = out1$res
res_global = out1$res_global
sample_data(ps)
# ancombc also supports importing data in phyloseq format
# tse_alt = agglomerateByRank(tse, "Family")
# pseq = makePhyloseqFromTreeSummarizedE

tab_lfc = res$lfc
head(tab_lfc)
#col_name = c("otu", "Intercept", "Fraction Total cells - Active ", "Fraction total DNA - Active", "Nodule-bulk soil", "Rhizo-bulksoil", 
#             "root-bulksoil")

col_name = c("otu", "Intercept", "Fraction Total cells - Active ", "Fraction total DNA - Active")

colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
tab_w$`Fraction Total cells - Active `


tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)
head(tab_p)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 7)
head(tab_q)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")
f<-filter(tab_diff, `Fraction Total cells - Active `==TRUE)

#####Format df from ANCOM######

#grab log fold change

tab_lfc<-tab_lfc %>% select(otu, `Fraction Total cells - Active `)
colnames(tab_lfc) <- c("otu", "lfc_total_active")
head(tab_lfc)

# grab pvalue
colnames(tab_p)
tab_p<-tab_p%>% select(otu, `Fraction Total cells - Active ` ) 
colnames(tab_p) <- c("otu","pvalue_total_active")
hist(tab_p$pvalue_total_active)
head(tab_p)

#grab adjusted p value
colnames(tab_q)
tab_q<-tab_q%>% select(otu, `Fraction Total cells - Active ` ) 
colnames(tab_q) <- c("otu","adj_pvalue_total_active")
hist(tab_q$adj_pvalue_total_active)

#grab TRUe FALSE
tab_diff<-tab_diff%>% select(otu, `Fraction Total cells - Active ` ) 
colnames(tab_diff) <- c("otu", "diff_total_active" )
tab_diff<-tab_diff%>% select(c(otu, diff_total_active))


#grab differentially abundant taxa
tab_diff<-tab_diff%>%filter(diff_total_active=="TRUE")
dim(tab_diff)

##otus table
row.names(df)
df$otu<-row.names(df)

# taxon table
taxon<-as.data.frame(tax_table(ps1))
df<-left_join(df, taxon)
df

#add ancom values
head(df)
df1<- df %>% select(-contains("R."))%>% 
  select(-contains("N."))%>% 
  select(-contains("E.")) %>%
  select(-contains("DNA_")) %>%

  left_join(., tab_p) %>%
  left_join(., tab_lfc) 

head(df1)

df1$pvalue_total_active



df_abund <- left_join(tab_diff, df2)
# negative lower in total than active 
# that confusing 
# multiplyu by -1 to make negative is lower in total
df2<-mutate(df2, lfc_total_active=lfc_total_active*-1)


# lfc verse the pvalue
ggplot(df1, aes(x=lfc_total_active , y=pvalue_total_active)) + 
  geom_jitter()+ #(stat="identity", position="identity") +
  #scale_fill_manual(values=c("#ab8af2" ,"#4c4b4d"))+
  #xlab("differentially abundant otus")+
  theme_bw()

colnames(df2)
df2$diff_total_active
# cut offs
df2<-df1 %>% filter(pvalue_total_active<.01, abs(lfc_total_active)>2)

# 146 differentially abundant taxa rhizosphere. 

df2

#plot 
#filter(df2,# abs(lfc_total_active)>3) 
ggplot(df2, aes(x=otu , y=lfc_total_active)) + 
  geom_bar(stat="identity", position="identity") +
    #scale_fill_manual(values=c("#ab8af2" ,"#4c4b4d"))+
  theme_bw()
# cute

#filter(df2, abs(lfc_total_active)>2) 
#filter(df2, df2$rhizo.total.sum>1)


filter(df2, df2$rhizo.total.mean>1) %>% filter(.,rhizo.bcat.mean>1)  %>%
ggplot(aes(x=otu , y=log10(rhizo.total.mean))) + 
  geom_bar(stat="identity", position="identity") +
  #xlab("abudance in active")+
  coord_flip()+
  theme_bw()
 

ggplot(df2, aes(x=rhizo.total.mean , y=lfc_total_active)) + 
  geom_jitter(stat="identity", position="identity") +
  xlab("abudance in total")+
  theme_bw()


# pvalue
df2$`Fraction Total cells - Active `

ggplot(df2, aes(x=rhizo.total.mean , y=log10(pvalue_total_active))) + 
  geom_jitter(stat="identity", position="identity") +
  xlab("abudance in total")+
  theme_bw()

df2$pvalue_total_active



# y = lfc totla active for all oty
# x = abundance
# color = differentailly abundant. 

df<-as.data.frame(t(otu_table(ps1)))
df$otu<-row.names(df)
# grab the the lfc and the pvalue and the true false
colnames(tab_diff) <- c("otu", "Intercept" ,"diff_total_active" ,"Fraction total DNA - Active")
tab_diff<-tab_diff%>% select(c(otu, diff_total_active))
tab_lfc
tab_q
tab<-left_join(tab_diff, tab_lfc) %>% left_join(., tab_q)
df2<-left_join(df, tab)

# sdd means
df2$rhizo.total.sum <-   rowSums(df2 %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df2$rhizo.total.mean <-   rowMeans(df2 %>% dplyr::select(contains("R.SYBR"))) %>% glimpse()
df2$rhizo.bcat.mean <-   rowMeans(df2 %>% dplyr::select(contains("R.POS"))) %>%   glimpse()
df2$rhizo.DNA.mean <-   rowMeans(df2 %>% dplyr::select(contains("R.DNA"))) %>%   glimpse()
colnames(df2)


# y = lfc totla active for all oty
# x = abundance
# color = differentailly abundant. 
dim(df2)
df2%>% filter(rhizo.total.mean<1000) %>% na.omit() %>%
ggplot( aes(x=rhizo.total.mean , y=lfc_total_active, color= diff_total_active)) + 
  geom_jitter(stat="identity", position="identity") +
  
  scale_color_manual(values=c("#f75477",
                             "#4c4b4d"))+
  xlab("abudance in total")+
  theme_bw()


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


#######regression active and total##############


filter(df) %>%
  ggplot(aes(x=rhizo.bcat.mean , y=rhizo.total.mean)) + 
  geom_jitter(stat="identity", position="identity") +
   theme_bw()

filter(df) %>%
  ggplot(aes(x=log10(rhizo.bcat.mean) , y=log10(rhizo.total.mean))) + 
  geom_jitter(stat="identity", position="identity") +
  geom_abline(slope=1)+
   theme_bw()

 df%>%
  ggplot(aes(x=nodule.xbcat.mean , y=nodule.total.mean)) + 
  geom_jitter(stat="identity", position="identity") +
   theme_bw()
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
 
 
 
 
