# flower Boncat
# March 2023
#Jennifer Harris

library(readxl)
library(tidyverse)
library(lubridate)
library(ggplot2)
library(forcats)
library(dplyr)
library(Hmisc)
library(lme4)
library(multcompView)
library(emmeans)
library(multcomp)
library(lmtest)

###-------oragnizing data Astrios----------

dir<-"C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/Flow cyto/Astrios"
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/Flow cyto/Astrios")
files<-list.files(path=dir, full.names = FALSE)

#import data
dat_csv <- lapply(files, read.csv)

#remove mean and SD rows
dat_csv<-lapply(dat_csv, function(k) subset(k, X!="Mean" & X!= "SD")) 
data1<-separate(data1, X, into = c(NA, NA, NA, NA, NA, "year", "month", "day", "sample"), extra= "merge")
dat_csv<-lapply(dat_csv, function(k) separate(k, X, into = c(NA, NA, NA, NA, NA, "year", "month", "day", "sample"), extra= "merge"))

# pick cols and rename
dat_csv<-lapply(dat_csv, function(k) subset(k, select = c(1:10)))
names=c("year", "month", "day", "sample", "total_count","singlet_count",
        "syto59_freq", "syto59_count", "BONCAT_freq", "BONCAT_count")
dat_csv <- lapply(dat_csv, function(x) setNames(x, names))

#turn into dataframe
df <- do.call("rbind", dat_csv)
#make dat column
df$Date<-as.Date(with(df,paste(year,month,day,sep="-")),"%Y-%m-%d")
write.csv(df, file= "combined.csv")
# added the metadat by hand in excel


###------import data --------
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/data")
#import data from fortessa and VYB flow cyto
df_fort<- read_excel("FlowCyto_Data_BONCAT_flower2021.xlsx")
counts<- read_excel("Flow cyto/cell_counts/Master_cell_counts.xlsx")
weights <- read_excel("Plant fitness/PlantFitness_Data_BONCAT_flower_Fall2020.xlsx", sheet = 2)
nodules<-read.csv("Plant fitness/nodule_surface_area.csv")
#import astrios data
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/Flow cyto/Astrios")
df_ast<-read.csv("combined_metadata_included.csv")
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
library(dplyr)
counts<-counts%>%
  mutate(cells_per_ul_diluted = green_pos_count / Volume_taken_ul)%>%
  mutate(cells_per_ul = cells_per_ul_diluted * dilution_1_XXX)%>%
  mutate(cells_per_gram_soil= ifelse( Fraction == "Bulk" | Fraction == "Rhizo", cells_per_ul*5000, NA))

# filtering for the cols we need from the surface area data
nodules<-nodules%>% select(c("Count", "Total.Area", "Plant", "Number", "Fraction"))
weights<-weights%>% mutate(Fraction = "Endo") %>% select("Plant", "Number", "abovegrnd_biomass_g", "Fraction")          

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
df<-df%>% mutate(Compartment=recode(Fraction, 'Bulk'='Bulk_Soil', 'Rhizo'='Rhizosphere','Endo'='Roots', 'Nod'='Nodule'))
df$Compartment<-factor(df$Compartment, levels = c("Bulk_Soil", "Rhizosphere", "Roots", "Nodule"))
df<-df%>% mutate(Plant=recode(Plant, 'A17'='Medicago', 'CLO'='Clover','PEA'='Pea', .default='soil'))

#-------set colors--------
mycols = c("#45924b", "#aac581", "#fffac9", "#f2a870", "#de425b")
mycols= c("#003f5c","#bc5090", "#ffa600")
#colors from 16s 
mycols1<- c("grey27", "#6eda0a")

#----set directory for figures------
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/data")

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

#Rhizo + bulk
svg(file="figures/rhizobulk_a17.svg",width = 4, height=4 )
rhizobulk%>% filter(Plant!="A17")%>%
  ggplot( aes(x=Fraction, y=cells_per_gram_soil)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.6,outlier.shape = NA , fill= mycols3[5])+
  scale_fill_manual(values = mycols3[5])+
  theme_minimal(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  scale_y_log10(limits= c(1E7, 1E9)) +
  #xlab("Plant")+
  ggtitle("Medicago")+
  ylab("Number Cells/ g soil")

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
# lil model of this. 
# I would consider these counts so poisson distribution seems appropriate
# um but poisson doesn;t fit very well becuase we have so many counts like 1 million counts is alot for her

hist(rhizobulk$cells_active_per_ul)
m1<-glm(rhizobulk$cells_active_per_ul~Plant+Fraction+Date+Plant*Fraction ,
        data=rhizobulk,
        family="poisson")
summary(m1)
plot(m1)

# um but poisson doesn;t fit very well becuase we have so many counts like 1 million counts is alot for her
# maybe i need to log transform

#log transformed#
#seems to be the best model for this
hist(log(rhizobulk$cells_active_per_ul))
# wooow so normal what a beauty
m1<-lm(log(rhizobulk$cells_active_per_ul)
  ~Plant+Fraction+Plant*Fraction, data = rhizobulk)
summary(m1)
anova(m1)
plot(m1)

# pretty strong interaction of fraction * plant
#breaking it up by plant for the analysis
#clo
hist(log(clo$cells_active_per_ul))
# pretty normal
clo<-rhizobulk%>%filter(Plant=="CLO")
m1<-lm(log(cells_active_per_ul)
       ~Fraction, data=clo)
summary(m1)
plot(m1)
# p = 1e-05
#a17
hist(log(a17$cells_active_per_ul))
# ehh
a17<-rhizobulk%>%filter(Plant=="A17")
m1<-lm(log(cells_active_per_ul)
            ~Fraction, data=clo)
summary(m1)
plot(m1)
# p is very tiny
#pea
hist(log(pea$cells_active_per_ul))
# ehh
pea<-rhizobulk%>%filter(Plant=="PEA")
m1<-lm(log(cells_active_per_ul)
       ~Fraction, data=pea)
summary(m1)
plot(m1)
# p is very tiny

# differences by plants
hist(log(rhizy$cells_active_per_ul))
# wooow so normal what a beauty
rhizy<-rhizobulk%>%filter(Fraction=="Rhizo")

m1<-lm(log(rhizy$cells_active_per_ul)
       ~Plant-1, data = rhizy)
summary(m1)
plot(m1)


#-------percent active figures-------

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/data")

#all plants
svg(file="figures/allfractions_percentwihtpea.svg",width = 13, height=5 )
df %>%
  filter(Plant!="soil", Dyes=="BONCAT-SYTO", flowcyto=="Fortessa")%>%
  ggplot( aes(x=Compartment, y=Percent_Boncat_pos)) +
  geom_boxplot(alpha=.7, outlier.shape = NA)+
  geom_jitter(width = .2)+
  scale_color_manual(values = mycols[c(1,2,3)])+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  facet_wrap(~Plant)+
  scale_y_log10()+
  xlab("compartment")+
  ylab("Percent Active")
  #ylim(0, 3E8)
dev.off()

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
hist(prop$n_Events_Cells)
y<-cbind(prop$n_Events_Boncat, prop$n_failures)

#mixed model
# it could possible that because each day I ran sampled on the flow Cyto they data is slightly different
# I should run a mixed model with allowing the intercep to vary for Date
# the intercept not the slope because I expect the relationship between the treatments to be the same but the baseline could vary.
# however, mixed models cannot have quasi binomial distribution. 
  m1<-glmer(y~Plant+Fraction+Plant*Fraction+(1|Date),
            data=prop,
            family="binomial")
  summary(m1)
  plot(m1)  
  ## residual deviance is higher than the degrees of freedom = model is over disposed
## these should be equal. so will will do a quasi binomial instead
# however, mixed models cannot have quasi binomial distribution. 
# I thought about including date as a fixed effect in my quasibinomial model, but date is very correlated with fraction. 
# quasibinomial glm

# glm on success and failures
  
  m4<-glm(data= prop, y~Plant, family = quasibinomial)
    #plot(m4)
    anova(m4, test= "LRT")
    m5<-glm(data= prop, y~Plant+Fraction, family = quasibinomial)
    #plot(m4)
    anova(m5, test= "LRT")
    lrtest(m5, m4)
    
##
  #             Df Deviance Resid. Df Resid. Dev        F    Pr(>F)    
  #   NULL                        82     379136                       
  #   Plant     4   366488        78      12648 406.4890 < 2.2e-16 ***
  #   Fraction  3     6564        75       6084   9.7079 1.729e-05 ***
  #    ---
  #     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
 summary(m4)
# Coefficients:
  #                 Estimate Std. Error t value Pr(>|t|)    
  #   PlantClover    -5.0073     0.9520  -5.260 1.32e-06 ***
  #   PlantMedicago  -6.4578     1.1370  -5.680 2.43e-07 ***
  #   PlantPea       -8.6241     1.3074  -6.596 5.25e-09 ***
  #   Plantsoil      -6.9075     5.6773  -1.217  0.22754    
  #   FractionRhizo  -0.3687     1.6752  -0.220  0.82639    
  #   FractionEndo    2.8519     0.9962   2.863  0.00544 ** 
  #   FractionNod     2.6623     0.9819   2.711  0.00831 ** 
 
 m4<-glm(data= prop, y~Plant+Fraction+Plant*Fraction-1, family = quasibinomial)
 #plot(m4)
 anova(m4, test= "LRT")
 summary(m4)
 
 m3 <- glm(data= prop, y~Plant+Fraction, family = binomial)
 lrtest(m3)
 glm
 anova(m3, test= "Chisq")
 
 
 
 
  
# just look at this for each plant
# need to find a good post hoctest
  clo<-prop%>% filter(Plant == "Clover")
  y<-cbind(clo$n_Events_Boncat, clo$n_failures)
  Fraction = clo$Fraction
  
  m3<-glm(y~Fraction, quasibinomial)
  #plot(m3)
  
  # do you just subset it all the way and run liley hood ratio test?
  anova(m3, test= "LRT")
  summary(m3)
  
  
  
  #pea
  pea<-prop%>% filter(Plant == "Pea")
  y<-cbind(pea$n_Events_Boncat, pea$n_failures)
  
  m3<-glm(y~pea$Fraction, quasibinomial)
  anova(m3, test= "F")
  summary(m3)
  plot(m3)
  
  #A17
  A17<-prop%>% filter(Plant == "Medicago")
  y<-cbind(A17$n_Events_Boncat, A17$n_failures)
  
  Fraction = A17$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m3, test= "F")
  summary(m2)
  plot(m2)
  

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
