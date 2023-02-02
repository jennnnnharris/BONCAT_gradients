# flower Boncat
# June 2021 
# project update
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
#data<-left_join(data, nodules, by= c("Plant", "Number", "Fraction"))
#data<-left_join(data, weights, by = c("Plant", "Number", "Fraction"))
#calculate number of active cells
ul_PBS = 10000

df<-df%>%mutate(cells_active_per_ul= Percent_Boncat_pos*cells_per_ul,
                    cells_active_g_soil = Percent_Boncat_pos*cells_per_gram_soil,
                    cells_per_plant= cells_per_ul*ul_PBS,
                    Prop_Active = Percent_Boncat_pos/100)

df$Fraction<-factor(df$Fraction, levels = c("Bulk", "Rhizo", "Endo", "Nod"))


#-------set colors--------
mycols = c("#45924b", "#aac581", "#fffac9", "#f2a870", "#de425b")

#-------n cells figures -------------

#bulk + rhizo by plant

 rhizobulk<- df %>%
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
  facet_wrap(vars(flowcyto))+
  xlab("Plant")+
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
# nodules + roots
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
#nod roots clover
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

###-------- linear model for n cells------
m1<-lm(cells_per_gram_soil~Plant+Fraction+Plant*Fraction, data = rhizobulk)
summary(m1)
plot(m1)

###-------log transformed-------#
#seems to be the best model for this
rhizobulk<-rhizobulk%>% filter(Plant!="PEA", Dyes=="BONCAT-SYTO")%>% 
  mutate(log_cell_g_soil = log(cells_per_gram_soil), Date = as.factor(Date))
m1<-lm(log_cell_g_soil~Plant+Fraction+Date+Plant*Fraction, data = rhizobulk)
summary(m1)
plot(m1)

###--------mixed model for n cells---------#####

cells_soil_mixed<- lmer(log(cells_per_gram_soil) ~ Fraction + Plant + (0|flowcyto), data = rhizobulk)
summary(cells_soil_mixed)
plot(cells_soil_mixed)
confint(cells_soil_mixed)
# varience of random effect is close to zero so this probably isn't the best model for the data.

#-----n active cells figures ---------------

#bulk
svg(file="figures/bulk_n_active_cells.avg",width = 6, height=6 )
df %>%
  filter(Dyes == "BONCAT-SYTO") %>%
  filter(Fraction == "Bulk") %>%
  ggplot( aes(x=Plant, y=n_Events_Boncat)) +
  geom_jitter(width = .2)+
  geom_boxplot(alpha=.5, fill = mycols[3], outlier.shape = NA)+
  #scale_fill_manual(values = mycols)+
  theme_bw(base_size = 20, )+
  #theme(axis.text.x = element_text(angle=60, hjust=1))+
  #facet_wrap(~Fraction)+
  xlab("Plant")+
  #ylab("Number Active Cells/ gram of soil")+
  #ylim(0,3E8)+
  #scale_y_log10() + 
  ggtitle("Bulk Soil")
dev.off()

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

#-------percent active figures-------

#all plants
svg(file="figures/allfractions_percent.svg",width = 6, height=6 )
data %>%
  filter(Plant!="Soil")%>%
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

#log scale
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

#clover
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

#a17
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

#-------binomial model for active cells ---------
# our data of boncat pos is a proprotion
# df frame with # failures # successes and proportion of each # i also got rid of an outlier past cook's distace
prop<-df %>%
  filter(Plant!="Soil", Dyes=="BONCAT-SYTO", ID!="beads") %>%
  #filter(Plant!="PEA")%>%
  mutate(n_failures = n_Events_Cells-n_Events_Boncat)
  #mutate(Active = ifelse(Prop_Active<.0001, 0, 1 ))
# data exploration
hist(prop$Prop_Active)
hist(prop$n_Events_Cells)

y<-cbind(prop$n_Events_Boncat, prop$n_failures)

#----------binomial model--------- 

m1<-glm(prop$Prop_Active~Plant+Fraction+Date,
          data=prop,
          family="binomial")
summary(m1)
plot(m1)  

m3<-glm(data= prop, Prop_Active~Plant+Fraction+Date-1, family=binomial, weights = n_Events_Cells)
summary(m3)

##------Poisson model-------

m1<-glm(prop$n_Events_Boncat~Plant+Fraction+Date+Plant*Fraction ,
        data=prop,
        family="poisson")
summary(m1)


#----------mixed model glmme-------
# it could possible that because each day I ran sampled on the flow Cyto they data is slightly different
# I should run a mixed model with allowing the intercep to vary for Date
# the intercept not the slope because I expect the relationship between the treatments to be the same but the baseline could vary.
# Random effect = Date
m1<-glmer(Prop_Active~Fraction +Plant+ Plant*Fraction + (1|Date),data = prop, binomial)
summary(m1)
plot(m1)
# modeling data as successes and failures

  y<-cbind(prop$n_Events_Boncat, prop$n_failures)
  m1<-glmer(y~Plant+Fraction+Plant*Fraction+(1|Date),
            data=prop,
            family="binomial")
  summary(m1)
  plot(m1)  
## residual deviance is higher than the degrees of freedom = model is over disposed
## these should be equal. so will will do a quasi binomial instead which u can't do in glmer
  
#------quasibinomial glm ----------

# on proportion data
  m2<-glm(Prop_Active~Plant+Fraction+Date, data=prop, quasibinomial)
  summary(m2)
  plot(m2)
  
  anova(m2, test= "F")  # significant interaction
  anova(m2, test= "Chisq")  

# weighted glm instead
  m3<-glm(data= prop, Prop_Active~Plant+Fraction+Date-1, family=quasibinomial, weights = n_Events_Cells)
  summary(m3)
  plot(m3)
  anova(m3, test= "F")

# glm on succeses and failures
  
  m4<-glm(data= prop, y~Plant+Fraction+Plant*Fraction+Date -1, family = quasibinomial)
  summary(m4)  
  plot(m4)
    
  # 2 step glm / hurdle model
  # glm for active verse non active

  
  Active<-filter(prop, Active=="1")
  m5<-glm(data=Active, Prop_Active~Plant*Fraction, family = binomial, weights = n_Events_Cells)
  summary(m5)
  
  #poisson
  m1<-glm(n_Events_Boncat~Fraction+ Plant + Plant*Fraction,  family= "poisson",
          offset(n_Events_Cells), data= prop)
  summary(m1)
  with(m1, cbind(res.deviance = deviance, df = df.residual,
                 p = pchisq(deviance, df.residual, lower.tail=FALSE)))
  
  # the P value i very low indicating that the data doesn't fit the model well.
  
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


#------some other models-------
