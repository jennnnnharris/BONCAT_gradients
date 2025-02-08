# Microbial activity analysis
# last updated 7 December 2024
# Jennifer Harris

library(readxl)
library(tidyverse)
library(dplyr)
library(lubridate)

###------import data --------

setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")
#import data from fortessa 
df<- read_excel("FlowCyto_data.xlsx")


#------cleaning data-------
#make date column date format
df$Date<-as.Date(df$Date , format = "%y%m$d")
df$Date<-ymd(df$Date)

#add dyes column
Dyes <-c()
Dyes <-ifelse(grepl("W",df$ID), 'unstained', 'BONCAT-SYTO')
df$Dyes <- Dyes
head(df)


### making levels a factor for plotting
df$Compartment<-factor(df$Compartment, levels = c("Bulk", "Rhizo", "Root", "Nod"))


#-------set colors--------
mycols= c("#003f5c","#bc5090", "#ffa600")

#----set directory for figures------
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Fig3_activity_diversity")

#-------percent active figures-------

svg(file="percent_active.svg",width = 5, height=4 )
df %>%
  filter(Treatment=="HPG")%>%
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
soil <- df%>%
  filter(Treatment == "HPG") %>%
  filter(Compartment!="Nod", Compartment!="Root") 

soil %>%  group_by(Compartment) %>%
  summarise(meanpercent=mean(Percent_Boncat_pos),
            SD= sd(Percent_Boncat_pos))

# avg % active in plant?
endosphere<-df%>%
  filter(Treatment == "HPG") %>%
  filter(Compartment!="Rhizo", Compartment!="Bulk") 

endosphere%>%
  group_by(Compartment) %>%
  summarise(meanpercent=mean(Percent_Boncat_pos),
            SD= sd(Percent_Boncat_pos))


#---------model for percent active cells ---------
# our data of BONCAT activity is a proportion
# to model this data some functions require a df frame with # failures # successes and proportion of each 

prop<-df %>%
  filter(Plant!="Soil", Treatment!="Methione",Dyes=="BONCAT-SYTO", ID!="beads", machine=='Fortessa') %>%
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
prop
# subset to compare compartment
  y<-cbind(prop$n_Events_Boncat, prop$n_failures)
  Fraction = prop$Compartment
  m2<-glm(y~prop$Compartment, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse rhizosphere
  prop1<-prop%>% filter(Compartment!="Endo") %>% filter(Compartment!="Nod")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Compartment = prop1$Compartment
  m2<-glm(y~Compartment, quasibinomial)
  anova(m2, test= "LRT")
    
# bulk verse endo
  prop1<-prop%>% filter(Compartment!="Rhizo") %>% filter(Compartment!="Nod")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Compartment = prop1$Compartment
  m2<-glm(y~Compartment, quasibinomial)
  anova(m2, test= "LRT")
  
# bulk verse nodules
  prop1<-prop%>% filter(Compartment!="Rhizo") %>% filter(Compartment!="Root")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Compartment = prop1$Compartment
    m2<-glm(y~Compartment, quasibinomial)
  anova(m2, test= "LRT")  

# rhizo verse endo
  prop1<-prop%>% filter(Compartment!="Bulk") %>% filter(Compartment!="Nod")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Compartment = prop1$Compartment
  m2<-glm(y~Compartment, quasibinomial)
  anova(m2, test= "LRT") 
  
# rhizo verse nodules
  prop1<-prop%>% filter(Compartment!="Bulk") %>% filter(Compartment!="Root")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Compartment = prop1$Compartment
  m2<-glm(y~Compartment, quasibinomial)
  anova(m2, test= "LRT") 
  
# endo verse nodules  
  prop1<-prop%>% filter(Compartment!="Bulk") %>% filter(Compartment!="Rhizo")
  y<-cbind(prop1$n_Events_Boncat, prop1$n_failures)
  Compartment = prop1$Compartment
  m2<-glm(y~Compartment, quasibinomial)
  anova(m2, test= "LRT") 
  

  
#-----cell counts-----

  
#import data

setwd("C:/Users/Jenn/OneDrive - The Pennsylvania State University/Documents/Github/BONCAT_gradients/data")
#import data from fortessa 
df<- read_excel("Cell_counts.xlsx")
  
#------data cleaning------
head(df)
### making levels a factor for plotting
df$Compartment<-factor(df$Compartment, levels = c("Bulk", "Rhizo", "Root", "Nodule"))
### remove unstained samples and flow cytometer controls. We will not be plotting these.
df<-df %>% filter(df$Stained=="1") 
# remove  flow cytometer controls. We will not be plotting these.
df<-df %>% filter(df$Plant!="CTL")
#make cell per ge soil numeric
df$approx_cells_per_g_soil<-as.numeric(df$approx_cells_per_g_soil)
#we are just going to focus on bulk soil and rhizo sphere
df<-df %>%filter(Compartment!= "Root" & Compartment!= "Nodule")

  
#set directory for figures
setwd("C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_gradients/Manuscript/figures/Supplement")
  
#-------cell counts figures-------
  
  svg(file="Cell_counts_per_g.svg",width = 5, height=4 )
  df %>%
    filter(Compartment!= "Root" & Compartment!= "Nodule")%>%
    ggplot(aes(x=Compartment, y=approx_cells_per_g_soil)) +
    geom_boxplot(alpha=.7, outlier.shape = NA)+
    geom_jitter(width = .1, size=1)+
    theme_classic(base_size = 14)+
    theme(axis.text.x = element_text(angle=60, hjust=1))+ 
    theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
          legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "none")+
    ylab("approximate cells / gram soil")
  dev.off()

  svg(file="Cell_counts.svg",width = 5, height=4 )
  df%>%
    filter(Compartment!= "Root" & Compartment!= "Nodule")%>%
    ggplot(aes(x=Compartment, y=cells_per_ul_sample)) +
    geom_boxplot(alpha=.7, outlier.shape = NA)+
    geom_jitter(width = .1, size=1)+
    theme_classic(base_size = 14)+
    theme(axis.text.x = element_text(angle=60, hjust=1))+ 
    theme(axis.text.x = element_text(angle=60, hjust=1, size = 14),
          legend.text = element_text(size = 14), axis.title.y =  element_text(size = 14),
          axis.text.y = element_text(size = 14),
          legend.position = "none")+
  ylab("cells / µL")
  dev.off()
  
  
#----summary statistics-----  
  
# avg % active in soil?
df %>%  group_by(Compartment) %>%
    summarise(mean_cell_count=mean(cells_per_ul_sample),
              SD_cell_count= sd(cells_per_ul_sample),
              mean_cell_per_g=mean(approx_cells_per_g_soil),
              SD_cell_per_g= sd(approx_cells_per_g_soil)
              )
  
#----anova----
  
m1<-lm(data=df, cells_per_ul_sample~Compartment )  
summary(m1)  
