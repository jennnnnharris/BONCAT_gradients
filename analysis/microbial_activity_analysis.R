# flower Boncat
# March 2023
#Jennifer Harris

#library(readxl)
#library(tidyverse)
#library(dplyr)
#library(lubridate)
#library(Hmisc) 
#library(lme4) # for generalized linear models
#library(multcompView)
#library(emmeans)
#library(multcomp)
#library(lmtest)


###------import data --------

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/")
#import data from fortessa and VYB flow cyto
df_fort<- read_excel("FlowCyto_Data_BONCAT_flower2021.xlsx")
counts<- read_excel("Flow cyto/cell_counts/Master_cell_counts.xlsx")
weights <- read_excel("Plant fitness/PlantFitness_Data_BONCAT_flower_Fall2020.xlsx", sheet = 2)
nodules<-read.csv("Plant fitness/nodule_surface_area.csv")


#import astrios data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Data/Flow cyto/")
df_ast<-read.csv("Astrios/combined_metadata_included.csv")
#------cleaning data-------
# remove rows without BONCAT dyes + controls
df_ast<-subset(df_ast, Dyes!= "SYBR- no click")
# add column for flow cyto machine
flowcyto<- rep("Astrios", length(df_ast$Date))
Plant<- rep("CLO", length(df_ast$Date))
df_ast<-cbind(df_ast, flowcyto)
df_ast<- cbind(df_ast, Plant)
#order columns
df_ast<-subset(df_ast, select=-c(year, month, day, singlet_count, total_count))
df_ast<-df_ast[,c(1,11,2,4,10,3, 6:9,5) ]
df_ast$Date<-mdy(df_ast$Date)  

#make date column date format
df_fort$Date<-as.Date(df_fort$Date , format = "%y%m$d")
df_fort$Date<-ymd(df_fort$Date)
#add flow cyto column
flowcyto<- rep("Fortessa", length(df_fort$Date))
df_fort<-cbind(df_fort, flowcyto)
Dyes <-c()
Dyes <-ifelse(grepl("W",df_fort$ID), 'unstained', 'BONCAT-SYTO')
Dyes[grepl("M",df_fort$ID)] = "meth"


df_fort<-cbind(df_fort, Dyes)
#remove unnesscary columns
df_fort<-subset(df_fort, select=-c(Project,Incubation,Number, Volume_filtered_ul, Spin_Speed, COUNT, Treatment, Run_number))
colnames(df_fort)
df_fort<-df_fort[,c(1:4,10,11, 5:9)]
#make column names the same
colnames(df_ast)<-colnames(df_fort)
# combine
df<-rbind(df_ast, df_fort)


#--------calculating cells/ ul and cells/ g soil-------
counts<-counts%>%
  mutate(cells_per_ul_diluted = green_pos_count / Volume_taken_ul)%>%
  mutate(cells_per_ul = cells_per_ul_diluted * dilution_1_XXX)%>%
  mutate(cells_per_gram_soil= ifelse( Fraction == "Bulk" | Fraction == "Rhizo", cells_per_ul*5000, NA))

    
#filter data, taking just the row and cols we need
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
df<-completeFun(df, 3)
counts<-counts[c(2,4,8:13)]

#merging boncat data and count data
df<-left_join(df, counts, by = c("ID", "Plant", "Fraction"))

#calculate number of active cells
ul_PBS = 10000

df<-df%>%mutate(cells_active_per_ul= Percent_Boncat_pos*cells_per_ul,
                    cells_active_g_soil = Percent_Boncat_pos*cells_per_gram_soil,
                    cells_per_plant= cells_per_ul*ul_PBS,
                    Prop_Active = Percent_Boncat_pos/100)

### recoding some names to be easier to understand
df$Fraction<-factor(df$Fraction, levels = c("Bulk", "Rhizo", "Endo", "Nod"))
df<-df%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
df$Compartment<-factor(df$Compartment, levels = c("Bulk Soil", "Rhizosphere", "Roots", "Nodule"))
df<-df%>% mutate(Plant=recode(Plant, 'A17'='Medicago', 'CLO'='Clover','PEA'='Pea', .default='soil'))
df$Fraction
df$Compartment

#-------set colors--------
mycols = c("#45924b", "#aac581", "#fffac9", "#f2a870", "#de425b")
mycols= c("#003f5c","#bc5090", "#ffa600")

#----set directory for figures------

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig3_activity_diversity")
#-------n cells figures -------------

#bulk + rhizo by plant

 rhizobulk<- df %>%
  filter(Plant!="Soil", Fraction!="Nod", Fraction!="Endo", !is.na(Fraction))
  rhizobulk$Fraction<-factor(rhizobulk$Fraction, levels = c("Bulk", "Rhizo"))
  rhizobulk<-subset(rhizobulk, Prop_Active<.2)
  rhizobulk<-subset(rhizobulk, cells_active_g_soil<1.5e+09)
  rhizobulk<-subset(rhizobulk, cells_active_g_soil>0)

svg(file="figures/bulkrhizo.svg",width = 4, height=4 )
  rhizobulk %>%
  ggplot(aes(x=Fraction, y=cells_per_gram_soil)) +
  geom_jitter(width = .2, size=1 )+
  geom_boxplot(alpha=.5, fill = "grey", outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(vars(flowcyto))+
  ylab("Number Cells/g soil")+
  scale_y_log10(limits= c(1E7, 1E9))
  
dev.off()

###--------  model for n cells------
#linear
m1<-lm(cells_per_gram_soil~Plant+Fraction+Plant*Fraction, data = rhizobulk)
summary(m1)
plot(m1)

#log
#seems to be the best model for this
rhizobulk<-rhizobulk%>% filter(Plant!="PEA", Dyes=="BONCAT-SYTO")%>% 
  mutate(log_cell_g_soil = log(cells_per_gram_soil), Date = as.factor(Date))
m1<-lm(log_cell_g_soil~Plant+Fraction+Date+Plant*Fraction, data = rhizobulk)
summary(m1)
plot(m1)

#mixed model
cells_soil_mixed<- lmer(log(cells_per_gram_soil) ~ Fraction + Plant + (0|flowcyto), data = rhizobulk)
summary(cells_soil_mixed)
plot(cells_soil_mixed)
confint(cells_soil_mixed)
# varience of random effect is close to zero so this probably isn't the best model for the data.

#-----n active cells figures ---------------

#rhizo bulk
svg(file="figures/bulk_n_active_cells.svg",width = 9, height=6 )
df%>%
  filter(Dyes == "BONCAT-SYTO", flowcyto=="Fortessa",Plant!="soil" ) %>%
  filter(Compartment!="Nodule", Compartment!="Roots") %>%
  ggplot( aes(x=Fraction, y=cells_active_g_soil)) +
  geom_jitter(width = .2, show.legend = FALSE)+
  geom_boxplot(alpha=.5, outlier.shape = NA, show.legend = FALSE)+
  #scale_fill_manual(values = mycols1)+
  theme_bw(base_size = 20, )+
  #theme(axis.text.x = element_text(angle=60, hjust=1))+
  facet_wrap(~Plant)+
  #xlab("Plant")+
  ylab("Number Active Cells/ gram of soil")+
  #ylim(0,3E8)+
  scale_y_log10()
  #ggtitle("Bulk Soil")
dev.off()

#--------model n active cells---------
#log transformed#
#seems to be the best model for this
hist(log(rhizobulk$cells_active_per_ul))
# wooow so normal what a beauty
m1<-lm(log(rhizobulk$cells_active_per_ul)
  ~Plant+Fraction+Plant*Fraction, data = rhizobulk)
anova(m1)
plot(m1)

# pretty strong interaction of fraction * plant
#breaking it up by plant for the analysis
#clo
hist(log(clo$cells_active_per_ul))
# pretty normal
clo<-rhizobulk%>%filter(Plant=="Clover")
m1<-lm(log(clo$cells_active_per_ul)
       ~clo$Fraction)
anova(m1)
# this code below does the same thing
a<-aov(log(clo$cells_active_per_ul)
    ~clo$Fraction)
summary(a)


#-------percent active figures-------
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig3_activity_diversity")

#clover
svg(file="percent_active.svg",width = 5, height=4 )
df %>%
  filter(Plant=="Clover", Dyes=="BONCAT-SYTO", flowcyto=="Fortessa")%>%
  ggplot( aes(x=Compartment, y=Percent_Boncat_pos)) +
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(width = .1, size=1)+
  scale_color_manual(values = mycols[c(1,2,3)])+
  theme_classic(base_size = 14)+
  theme(axis.text.x = element_text(angle=60, hjust=1))+ 
  theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
        legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
          axis.text.y = element_text(size = 14),
            legend.position = "none")+
  scale_y_log10()+
  xlab("")+
  ylab("Percent Active cells")
#ylim(0, 3E8)
dev.off()

# avg % active in soil?
df%>%
  filter(Dyes == "BONCAT-SYTO", flowcyto=="Fortessa", Plant=="Clover") %>%
  filter(Compartment!="Nodule", Compartment!="Roots") %>% #group_by(Compartment) %>%
  summarise(meanpercent=mean(Percent_Boncat_pos),
            SD= sd(Percent_Boncat_pos))

df%>%
  filter(Dyes == "BONCAT-SYTO", flowcyto=="Fortessa", Plant=="Clover") %>%
  filter(Compartment!="Rhizosphere", Compartment!="Bulk Soil") %>% #group_by(Compartment) %>%
  summarise(meanpercent=mean(Percent_Boncat_pos),
            SD= sd(Percent_Boncat_pos))

#---------model for percent active cells ---------
# our data of boncat activity is a proportion
# to model this data some functions require a df frame with # failures # successes and proportion of each 
# i also got rid of an outlier past cook's distace
prop<-df %>%
  filter(Plant!="Soil", Dyes=="BONCAT-SYTO", ID!="beads", flowcyto=='Fortessa') %>%
  mutate(n_failures = n_Events_Cells-n_Events_Boncat)
  #mutate(Active = ifelse(Prop_Active<.0001, 0, 1 ))
# data exploration
hist(prop$Prop_Active)
hist(log(prop$Prop_Active))
hist(prop$n_Events_Cells)
y<-cbind(prop$n_Events_Boncat, prop$n_failures)

# mixed models cannot have quasi binomial distribution. 
# I need a quasibinomial model because residual deviance is higher than the degrees of freedom = model is over disposed
# I thought about including date as a fixed effect in my quasibinomial model, but date is very correlated with fraction. 
# quasibinomial glm on success and failures

 m1<-glm(data= prop, y~Plant+Fraction+Plant*Fraction, family = quasibinomial)
 m1
 
 anova(m1, test= "LRT")  
 anova(m1, test= "Chisq")

# Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
# NULL                              81    20968.2              
# Plant           3   8320.2        78    12648.0 < 2.2e-16 ***
# Fraction        3   6564.4        75     6083.6 1.177e-15 ***
# Plant:Fraction  6   1283.1        69     4800.5   0.02756 * 
 
# significant interaction 
# break into each plant 
# clover  
clo<-prop%>% filter(Plant == "Clover")
  y<-cbind(clo$n_Events_Boncat, clo$n_failures)
  Fraction = clo$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse rhizosphere
  clo1<-clo%>% filter(Fraction!="Endo") %>% filter(Fraction!="Nod")
  y<-cbind(clo1$n_Events_Boncat, clo1$n_failures)
  Fraction = clo1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
    
# bulk verse endo
  clo1<-clo%>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Nod")
  y<-cbind(clo1$n_Events_Boncat, clo1$n_failures)
  Fraction = clo1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse nodules
  clo1<-clo%>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Endo")
  y<-cbind(clo1$n_Events_Boncat, clo1$n_failures)
  Fraction = clo1$Fraction
    m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")  

# rhizo verse endo
  clo1<-clo%>% filter(Fraction!="Bulk") %>% filter(Fraction!="Nod")
  y<-cbind(clo1$n_Events_Boncat, clo1$n_failures)
  Fraction = clo1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
# rhizo verse nodules
  clo1<-clo%>% filter(Fraction!="Bulk") %>% filter(Fraction!="Endo")
  y<-cbind(clo1$n_Events_Boncat, clo1$n_failures)
  Fraction = clo1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
# endo verse nodules  
  clo1<-clo%>% filter(Fraction!="Bulk") %>% filter(Fraction!="Rhizo")
  y<-cbind(clo1$n_Events_Boncat, clo1$n_failures)
  Fraction = clo1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
######pea
  pea<-prop%>% filter(Plant == "Pea")
  y<-cbind(pea$n_Events_Boncat, pea$n_failures)
  Fraction = pea$Fraction

  # bulk verse rhizosphere
  pea1<-pea %>% filter(Fraction!="Endo") %>% filter(Fraction!="Nod")
  y<-cbind(pea1$n_Events_Boncat, pea1$n_failures)
  Fraction = pea1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
  # bulk verse endo
  pea1<-pea %>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Nod")
  y<-cbind(pea1$n_Events_Boncat, pea1$n_failures)
  Fraction = pea1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
  # bulk verse nodules
  pea1<-pea %>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Endo")
  y<-cbind(pea1$n_Events_Boncat, pea1$n_failures)
  Fraction = pea1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")  
  
  # rhizo verse endo
  pea1<-pea %>% filter(Fraction!="Bulk") %>% filter(Fraction!="Nod")
  y<-cbind(pea1$n_Events_Boncat, pea1$n_failures)
  Fraction = pea1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
  # rhizo verse nodules
  pea1<- pea %>% filter(Fraction!="Bulk") %>% filter(Fraction!="Endo")
  y<-cbind(pea1$n_Events_Boncat, pea1$n_failures)
  Fraction = pea1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
  # endo verse nodules  
  pea1<-pea %>% filter(Fraction!="Bulk") %>% filter(Fraction!="Rhizo")
  y<-cbind(pea1$n_Events_Boncat, pea1$n_failures)
  Fraction = pea1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")   
  
  
#######A17
  A17<-prop%>% filter(Plant == "Medicago")
  y<-cbind(A17$n_Events_Boncat, A17$n_failures)
  Fraction = A17$Fraction
  
# bulk verse rhizosphere
  A171<-A17 %>% filter(Fraction!="Endo") %>% filter(Fraction!="Nod")
  y<-cbind(A171$n_Events_Boncat, A171$n_failures)
  Fraction = A171$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse endo
  A171<-A17 %>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Nod")
  y<-cbind(A171$n_Events_Boncat, A171$n_failures)
  Fraction = A171$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
  # bulk verse nodules
  A171<-A17 %>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Endo")
  y<-cbind(A171$n_Events_Boncat, A171$n_failures)
  Fraction = A171$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")  
  
  # rhizo verse endo
  A171<-A17 %>% filter(Fraction!="Bulk") %>% filter(Fraction!="Nod")
  y<-cbind(A171$n_Events_Boncat, A171$n_failures)
  Fraction = A171$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
  # rhizo verse nodules
  A171<- A17 %>% filter(Fraction!="Bulk") %>% filter(Fraction!="Endo")
  y<-cbind(A171$n_Events_Boncat, A171$n_failures)
  Fraction = A171$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
  # endo verse nodules  
  A171<-A17 %>% filter(Fraction!="Bulk") %>% filter(Fraction!="Rhizo")
  y<-cbind(A171$n_Events_Boncat, A171$n_failures)
  Fraction = A171$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 

#------Testing correlation of variables-------
  
data<-df[,c(1,3,6)]
data<-filter(data, Dyes=="BONCAT-SYTO")
data<-data[,c(1,2)]

binary<-data%>% mutate(binary=recode(Fraction, 'Endo'='plant', 'Nod'='plant',.default='soil'))
binary<-binary[-143, c(1,3)]

data<-data.table::as.data.table(data)
data<-data.table::dcast(data, Date~Fraction,fun.aggregate=length)
data<-data[,3:6] 
data<-as.matrix(data)
soil<-rowSums(data[-1,c(1,2)])
plant<-rowSums(data[-1,c(3,4)])

data2<-cbind(soil, plant)
# chi square test can be done if one of the categorical vars is binary (plant, soil)

test<-chisq.test(data2)
# X-squared = 54.12, df = 15, p-value = 2.509e-06
# p is less than .05 so variable are independent.
# we can normalize values on a line from 0 to 1 because chi square can because
# Pearson’s χ2 maximum value depends on the sample size and the size of the contingency table

library('DescTools')

ContCoef(data2, correct = FALSE)
ContCoef(data2, correct = TRUE)

library('rcompanion')

cramerV(data2)
cramerV(data2, bias.correct = TRUE)

install.packages("vcd")
library('vcd')

test<-assocstats(xtabs(~binary$Date + binary$binary))
test
#fishers exact test to see if date corr with sample type.

fisher.test(data2)



chisq.test(data)
