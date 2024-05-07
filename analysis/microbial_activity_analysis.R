# Microbial activity analysis
# May 2024
#Jennifer Harris

library(readxl)
library(tidyverse)
library(dplyr)
library(lubridate)

###------import data --------

setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")
#import data from fortessa and VYB flow cyto
df<- read_excel("FlowCyto_data.xlsx")


#------cleaning data-------
#make date column date format
df$Date<-as.Date(df_fort$Date , format = "%y%m$d")
df$Date<-ymd(df_fort$Date)

#add dyes column
Dyes <-c()
Dyes <-ifelse(grepl("W",df_fort$ID), 'unstained', 'BONCAT-SYTO')
df$Dyes <- Dyes
head(df)


### making levels a factor for plotting
df$Compartment<-factor(df$Compartment, levels = c("Bulk", "Rhizo", "Endo", "Nod"))


#-------set colors--------
mycols= c("#003f5c","#bc5090", "#ffa600")

#----set directory for figures------
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig3_activity_diversity")

#-------percent active figures-------

svg(file="percent_active.svg",width = 5, height=4 )
df %>%
  filter(df$Dyes=="BONCAT-SYTO")%>%
  ggplot(aes(x=Compartment, y=Percent_Boncat_pos)) +
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
  filter(Dyes == "BONCAT-SYTO") %>%
  filter(Compartment!="Nodule", Compartment!="Roots") %>% #group_by(Compartment) %>%
  summarise(meanpercent=mean(Percent_Boncat_pos),
            SD= sd(Percent_Boncat_pos))

# avg % active in plant?
df%>%
  filter(Dyes == "BONCAT-SYTO") %>%
  filter(Compartment!="Rhizosphere", Compartment!="Bulk Soil") %>% #group_by(Compartment) %>%
  summarise(meanpercent=mean(Percent_Boncat_pos),
            SD= sd(Percent_Boncat_pos))

#---------model for percent active cells ---------
# our data of BONCAT activity is a proportion
# to model this data some functions require a df frame with # failures # successes and proportion of each 

prop<-df %>%
  filter(Plant!="Soil", Dyes=="BONCAT-SYTO", ID!="beads", flowcyto=='Fortessa') %>%
  mutate(n_failures = n_Events_Cells-n_Events_Boncat)
  
y<-cbind(prop$n_Events_Boncat, prop$n_failures)

# mixed models cannot have quasi binomial distribution. 
# I need a quasibinomial model because residual deviance is higher than the degrees of freedom = model is over disposed
# I thought about including date as a fixed effect in my quasibinomial model, but date is very correlated with fraction. 
# quasibinomial glm on success and failures

 m1<-glm(data= prop, y~Compartment, family = quasibinomial)
 m1
 
 anova(m1, test= "LRT")  
 anova(m1, test= "Chisq")

# subset to compare compartment
  y<-cbind(prop$n_Events_Boncat, prop$n_failures)
  Fraction = prop$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse rhizosphere
  prop1<-prop%>% filter(Fraction!="Endo") %>% filter(Fraction!="Nod")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Fraction = prop1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
    
# bulk verse endo
  prop1<-prop%>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Nod")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Fraction = prop1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse nodules
  prop1<-prop%>% filter(Fraction!="Rhizo") %>% filter(Fraction!="Endo")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Fraction = prop1$Fraction
    m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT")  

# rhizo verse endo
  prop1<-prop%>% filter(Fraction!="Bulk") %>% filter(Fraction!="Nod")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Fraction = prop1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
# rhizo verse nodules
  prop1<-prop%>% filter(Fraction!="Bulk") %>% filter(Fraction!="Endo")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Fraction = prop1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
# endo verse nodules  
  prop1<-prop%>% filter(Fraction!="Bulk") %>% filter(Fraction!="Rhizo")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Fraction = prop1$Fraction
  m2<-glm(y~Fraction, quasibinomial)
  anova(m2, test= "LRT") 
  
