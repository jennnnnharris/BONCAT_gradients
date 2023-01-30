# flower Boncat
# June 2021 
# project update
#Jennifer Harris

library(readxl)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(forcats)
library(dplyr)
library(Hmisc)
library(lme4)
library(multcompView)
library(emmeans)
library(multcomp)



setwd("C:/Users/harri/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/data")
#import data
data<- read_excel("FlowCyto_Data_BONCAT_flower2021.xlsx")
counts<- read_excel("Flow cyto/cell_counts/Master_cell_counts.xlsx")
weights <- read_excel("Plant fitness/PlantFitness_Data_BONCAT_flower_Fall2020.xlsx", sheet = 2)
nodules<-read.csv("Plant fitness/nodule_surface_area.csv")
# calculating cells/ ul and cells/ g soil

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
data<-completeFun(data, 3)
counts<-counts[c(2,4,8:13)]

#merging boncat data and count data
data<-left_join(data, counts, by = c("ID", "Plant", "Fraction"))
#data<-left_join(data, nodules, by= c("Plant", "Number", "Fraction"))
#data<-left_join(data, weights, by = c("Plant", "Number", "Fraction"))
#calculate number of active cells
ul_PBS = 10000

data<-data%>%mutate(cells_active_per_ul= Percent_Boncat_pos*cells_per_ul,
                    cells_active_g_soil = Percent_Boncat_pos*cells_per_gram_soil,
                    cells_per_plant= cells_per_ul*ul_PBS,
                    Prop_Active = Percent_Boncat_pos/100)


#mycols = c("#003f5c", "#7a5195", "#ef5675", "#ffa600")
#mycols2 = c("#003f5c",   "#bc5090",   "#ffa600")

mycols = c("#45924b", "#aac581", "#fffac9", "#f2a870", "#de425b")


data$Fraction<-factor(data$Fraction, levels = c("Bulk", "Rhizo", "Endo", "Nod"))

hist(data$cells_per_ul, breaks = 30)
hist(data$cells_active_g_soil, breaks = 50)

###----- import data --------

dir<-"C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/Flow cyto/Astrios"

setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT/Data/Flow cyto/Astrios")
files<-list.files(path=dir, full.names = FALSE)

#import data
dat_csv <- lapply(files, read.csv)
data1<- read.csv(files[1])
data2<- read.csv(files[2])

#remove mean and SD rows
data1<-data1[!(data1$X=="Mean" | data1$X=="SD"),]
data2<-data2[!(data2$X=="Mean" | data2$X=="SD"),]
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


#-------n cells-------------

#bulk + rhizo by plant

 rhizobulk<- data %>%
  filter(Plant!="Soil", Fraction!="Nod", Fraction!="Endo", !is.na(Fraction))
rhizobulk$Fraction<-factor(rhizobulk$Fraction, levels = c("Bulk", "Rhizo"))  

svg(file="figures/bulkrhizo_clo.svg",width = 4, height=4 )
  rhizobulk %>% filter(Plant=="CLO")%>%
  ggplot(aes(x=Fraction, y=cells_per_gram_soil)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, fill = mycols[5], outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_minimal(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  facet_wrap(vars(Plant))+
  xlab("Plant")+
  ylab("Number Cells/g soil")+
  scale_y_log10(limits= c(1E7, 1E9))
  
dev.off()

rhizobulk<-rhizobulk%>% filter(Plant!="PEA")%>%
mutate(log_cell_g_soil = log(cells_per_gram_soil))
m1<-lm(log_cell_g_soil~Plant+Fraction+Plant*Fraction, data = rhizobulk)
summary(m1)
plot(m1)

a1<-aov(lm1$cells_per_gram_soil~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)


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

dev.off()

#t  test
lm1 <- data %>%
  filter(Plant!="Soil", Fraction!="Endo", Fraction!="Nod", Incubation=="H")
  t.test(cells_per_gram_soil~Fraction, data = lm1)
  

a1<-aov(lm1$cells_per_gram_soil~lm1$Fraction)
TukeyHSD(a1, 'lm1$Fraction', conf.level = .95)


nodendo<-data %>%
  filter(Plant!="Soil", Fraction=="Nod"|Fraction == "Endo")
  nodendo$Fraction<-factor(nodendo$Fraction, levels = c("Endo", "Nod"))
svg(file="figures/nodendo_a17.svg",width = 4, height=4 )
  nodendo%>% filter(Plant=="A17")%>%
  ggplot( aes(x=Fraction, y=cells_per_plant)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.6, fill = mycols[5], outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_minimal(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  facet_wrap(vars(Plant))+
  xlab("Plant")+
  ylab("Number Cells/Plant")+
  scale_y_log10(limits= c(1E6, 1E10))

dev.off()

svg(file="figures/nodendos_clo.svg",width = 4, height=4 )
nodendo%>% filter(Plant=="CLO")%>%
  ggplot( aes(x=Fraction, y=cells_per_plant)) +
  geom_jitter(width = .2, size=2 )+
  geom_boxplot(alpha=.5, fill = mycols3[5], outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_minimal(base_size = 22, )+
  ggtitle("clo")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(vars(Plant))+
  #xlab("Plant")+
  ylab("Number Cells/Plant")+
  scale_y_log10(limits= c(1E6, 1E10))
dev.off()


lm1 <- data %>%
  filter(Plant!="Soil", Fraction!="Endo", Fraction!="Nod", Incubation=="H", Plant!="PEA")
#t.test(cells_per_gram_soil~Plant, data = lm1)

a1<-aov(lm1$cells_per_gram_soil~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)


#t  test
lm1 <- data %>%
  filter(Plant!="Soil", Fraction!="Bulk", Fraction!="Rhizo", Incubation=="H")
t.test(cells_per_ul~Fraction, data = lm1)



#-----n active cells---------------

#bulk

svg(file="figures/bulk_n_active_cells.avg",width = 6, height=6 )
data %>%
  filter(Plant!="Soil", Fraction==("Bulk"), Treatment=="HPG")%>%
  ggplot( aes(x=Plant, y=cells_active_g_soil)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, fill = mycols[3], outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  xlab("Plant")+
  ylab("Number Active Cells/ gram of soil")+
  ylim(0,3E8)+
  #scale_y_log10() + 
  ggtitle("Bulk Soil")
dev.off()

lm1 <- data %>%
  filter(Plant!="Soil", Fraction=="Bulk", Incubation=="H")
#t.test(cells_active_g_soil~Plant, data = lm1)

a1<-aov(lm1$cells_active_g_soil~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)

#Rhizo
pdf(file="figures/rhizo_n_active_cells.pdf",width = 6, height=6 )
data %>%
  filter(Plant!="Soil", Fraction==("Rhizo"), Treatment=="HPG")%>%
  ggplot( aes(x=Plant, y=cells_active_g_soil)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, fill = mycols2,outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  xlab("Plant")+
  ylab("Number Active Cells/ gram soil")+
  ylim(0, 3E8)+
  ggtitle("Rhizosphere Soil")
dev.off()

lm1 <- data %>%
  filter(Plant!="Soil", Fraction=="Rhizo", Incubation=="H")
#t.test(cells_active_g_soil~Plant, data = lm1)

a1<-aov(lm1$cells_active_g_soil~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)


lm1 <- data %>%
  filter(Plant!="Soil", Fraction!="Endo", Fraction!="Nod", Incubation=="H")
#t.test(cells_per_gram_soil~Plant, data = lm1)

a1<-aov(lm1$cells_per_gram_soil~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)


#Rhizo + bulk
pdf(file="figures/rhizoandbulkactive.pdf",width = 6, height=6 )
data %>%
  filter(Plant!="Soil", Plant=="CLO", Fraction!="Endo", Fraction!="Nod", Incubation=="H")%>%
  ggplot( aes(x=Fraction, y=cells_active_g_soil)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  xlab("Plant")+
  ylab("Cells Active/ g soil")+
  scale_y_log10()
  ylim(0, 3E8)

dev.off()

#Rhizo + bulk
svg(file="figures/rhizoandbulkactive.svg",width = 3, height=3 )
data %>%
  filter(Plant!="Soil", Plant=="A17", Fraction!="Endo", Fraction!="Nod", Incubation=="H")%>%
  ggplot( aes(x=Fraction, y=cells_active_g_soil)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  xlab("Plant")+
  ylab("Cells Active/ g soil")+
  scale_y_log10()
ylim(0, 3E8)

dev.off()






#t  test
lm1 <- data %>%
  filter(Plant!="Soil", Fraction!="Endo", Fraction!="Nod", Incubation=="H")
t.test(cells_active_g_soil~Fraction, data = lm1)



#nod
pdf(file="figures/nod_n_active_cells_logscale.pdf",width = 6, height=6 )
data %>%
  filter(Plant!="Soil", Fraction=("Nod"), Treatment=="HPG")%>%
  ggplot( aes(x=Plant, y=cells_active_per_ul)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, fill = mycols2, outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  xlab("Plant")+
  ylab("N Active Cells/ Surface Area cm^2")+
  scale_y_log10()+
  #ylim(0, 4E5)+
  ggtitle("Nodule")
dev.off()
  

lm1 <- data %>%
  filter(Plant!="Soil", Fraction=="Nod", Incubation=="H")

a1<-aov(lm1$cells_active_per_ul~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)


#endo
pdf(file="figures/endo_nactive_cells.pdf",width = 6, height=6 )
data %>%
  filter(Plant!="Soil", Fraction==("Endo"), Treatment=="HPG")%>%
  ggplot( aes(x=Plant, y=cells_active_g_plant)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, fill = mycols2, outlier.shape = NA)+
  scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  xlab("Plant")+
  ylab("Number Active Cells/g plant")+
  #ylim(0, 4E5)+
  ggtitle("Endo")
dev.off()


lm1 <- data %>%
  filter(Plant!="Soil", Fraction=="Endo", Incubation=="H")

a1<-aov(lm1$cells_active_per_ul~lm1$Plant)
TukeyHSD(a1, 'lm1$Plant', conf.level = .95)

#is there an interaction with plant type?

# log transformed linear model 
colnames(data)

NE<-data%>%filter(Plant!="PEA", Plant!="Soil", Treatment=="HPG", Fraction!="Endo", Fraction!="Nod")%>%
  mutate(lg_cell_act_soil= log(cells_active_g_soil))
  
m1<-lm(lg_cell_act_soil~ Fraction + Plant+ Plant*Fraction, data = NE)
anova(m1)



NE<-data%>%filter(Plant!="PEA", Plant!="Soil", Treatment=="HPG", Fraction!="Endo", Fraction!="Nod")%>%
  mutate(lg_cell_act_soil= log(cells_active_g_soil))

m1<-lm(lg_cell_act_soil~ Fraction + Plant+ Plant*Fraction, data = NE)

a1<-aov(m1)
TukeyHSD(a1, conf.level = .95)


#-------percent-------

#all plants

svg(file="figures/allfractions_percent.svg",width = 6, height=6 )
data %>%
  filter(Plant!="Soil", Incubation=="H")%>%
  ggplot( aes(x=Fraction, y=Percent_Boncat_pos)) +
  geom_boxplot(alpha=.7, outlier.shape = NA, fill= mycols[4])+
  geom_jitter(width = .2)+
  theme_bw(base_size = 20, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  scale_y_log10()+
  xlab("Fraction")+
  ylab("Percent Active")
  #ylim(0, 3E8)
dev.off()



svg(file="figures/percentBONCATlog.svg",width = 14, height=6 )
data %>%
  filter(Treatment == "HPG", Plant!= "Soil", Plant!="PEA")%>%
  ggplot( aes(x=Fraction, y=Percent_Boncat_pos)) +
  geom_boxplot( alpha=.7, fill= mycols[4], outlier.shape = NA)+
  geom_jitter(size = 1.5, width = .2)+
  theme_minimal(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylim(0, 50)+
  xlab("Fraction")+
  ylab("% Boncat Active")+
  scale_y_log10(limits= c())+
  facet_wrap(vars(Plant))+
  ggtitle("Percent BONCAT active")
dev.off()


svg(file="figures/percentBONCATlog_clo.svg",width = 8, height=6 )
data %>%
  filter(Treatment == "HPG", Plant=="CLO")%>%
  ggplot( aes(x=Fraction, y=Percent_Boncat_pos)) +
  geom_boxplot( alpha=.7, fill= mycols[4], outlier.shape = NA)+
  geom_jitter(size = 1.5, width = .2)+
  theme_minimal(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylim(0, 50)+
  xlab("Fraction")+
  ylab("% Boncat Active")+
  scale_y_log10(limits= c())+
  facet_wrap(vars(Plant))
  #ggtitle("Percent BONCAT active")
dev.off()
#combined figure

svg(file="figures/percentBONCAT_log_a17.svg",width = 8, height=6 )
data %>%
  filter(Treatment == "HPG", Plant=="A17")%>%
  ggplot( aes(x=Fraction, y=Percent_Boncat_pos)) +
  geom_boxplot( alpha=.7, fill= mycols[4], outlier.shape = NA)+
  geom_jitter(size = 1.5, width = .2)+
  theme_minimal(base_size = 22, )+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  #ylim(0, 50)+
  xlab("Fraction")+
  ylab("% Boncat Active")+
  facet_wrap(vars(Plant))+
  scale_y_log10(limits= c())
dev.off()

#-------stats---------
# make a smaller df for tests
prop <- data %>%
  mutate(Plant_rep = paste0(Plant, Number))%>%
  filter(Plant!="Soil", Incubation=="H", Plant!="PEA")%>%
  select(Plant, Plant_rep, Fraction, Prop_Active,  n_Events_Cells, n_Events_Boncat)%>%
  mutate(n_failures = n_Events_Cells-n_Events_Boncat)%>%
  mutate(Active = ifelse(Prop_Active<.0001, 0, 1 ))
 

# mixed model
# Random effect = plantrep
# this would be a good way to account for if some plants are just more active than others...
# doesn't work super well here tho because then i basically get a dif random effect for each sample
# which doesn't make sense. 

lmer(Prop_Active~Fraction+ (Fraction|Plant_rep), data = prop)


# weighted glm instead

m3<-glm(data= prop, Prop_Active~Plant+Fraction+Plant*Fraction -1, family=quasibinomial, weights = n_Events_Cells)
summary(m3) 
plot(m3)

# 2 step glm / hurdle model

# glm for active verse non active
  m4<-glm(data= prop, Active~Plant+Fraction+Plant*Fraction -1, family = quasibinomial)
  summary(m4)

  Active<-filter(prop, Active=="1")
  m5<-glm(data=Active, Prop_Active~Plant*Fraction, family = binomial, weights = n_Events_Cells)
  summary(m5)
 
# random effect for plant  
  
glmm  
   
# vector of failures and successes
  y<-cbind(prop$n_Events_Boncat, prop$n_failures)
  
  m1<-glm(y~prop$Plant*prop$Fraction - 1, binomial)
  summary(m1)
  plot(m1)
  ## residual deviance is 3931 with 49 degrees of freedom.
  ## these should be equal. so will will do a quasi bonimial instead
  
  m1<-glm(y~prop$Fraction -1, quasibinomial)
  summary(m1)
  plot(m1)
  
  m2<-glm(y~prop$Plant+prop$Fraction+prop$Plant*prop$Fraction -1, quasibinomial)
  summary(m2)
  plot(m2)
  
  anova(m2, test= "F")  # significant interaction
  anova(m2, test= "Chisq")  


  
# just look at this for each plant
  clo<-prop%>% filter(Plant == "CLO")
  y<-cbind(clo$n_Events_Boncat, clo$n_failures)
  
  Fraction = clo$Fraction
  m3<-glm(y~Fraction, quasibinomial)
  summary(m3)
  plot(m3)
  
  anova(m3, test= "LR")
  
 
  marginal = emmeans(m3, ~Fraction)
  pairs(marginal, adjust = "tukey")

  
  summary(glht(m3,mcp(Fraction = "Tukey")))
  summary(glht(age, mcp(age="Tukey")))
  
  
  clo$Fraction<-factor(clo$Fraction, levels = c("Endo", "Nod", "Rhizo", "Bulk"))
  #pea
  pea<-prop%>% filter(Plant == "PEA")
  y<-cbind(pea$n_Events_Boncat, pea$n_failures)
  
  m2<-glm(y~pea$Fraction, quasibinomial)
  summary(m2)
  plot(m2)
  
  #A17
  A17<-prop%>% filter(Plant == "A17")
  y<-cbind(A17$n_Events_Boncat, A17$n_failures)
  
  Fraction = A17$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  summary(m2)
  plot(m2)
  
  anova(m2, test= "LR")
  
  summary(glht(m2,mcp(Fraction = "Tukey")))
  


#poisson
m1<-glm(Percent_Boncat_pos~Fraction+ Plant + Plant*Fraction,  family= "poisson", data= percent)
summary(m1)
with(m1, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))

# the P value i very low indicating that the data doesn't fit the model well.


