library(ggplot2)
library(dplyr)
library(vegan)
library(lme4)

#Read in dataset
df <- read.csv("./02_outdata/revisedSep18-25_backswimmer_development_dispersal.csv")


# Temperature on Dispersal ------------------------------------------------
#Exploring the data
disp <- subset(df, Disp_Count >= 0) #only 30 individuals survived to the dispersal assay

#Reviewing the data
  #table of individuals in dispersal propensity assay by temperature treatment
  disp_count <- table(disp$Treatment, disp$Disp_Count) #only 3 dispersal attempts
  
  #table of individuals that dispersed in dispersal ability assay by temperature
  disp_final <- table(disp$Treatment, disp$Dispersed) 
  sum(disp$Dispersed) #only 1 individual actually dispersed away from the experiment

#Takeaway - dispersal was extremely low in this experiment and the sample size is very small
  #Any analysis lacks power
 
#Generalized linear mixed effects model for temperature on dispersal attempts count
  #This uses a poisson distribution since the response is count data
dispcount_model <- glmer(Disp_Count ~ Treatment + (1|Block/Mesocosm), data = disp, family = poisson)
summary(dispcount_model)
  #Temperature has no significant effect on number of dispersal attempts

#Generalized linear mixed effects model for temperature on dispersal (absolute)
  #This uses binomial distribution since the response is binary (1 = dispersed, 0 = stayed) 
disp_model <- glmer(Dispersed ~ Treatment + (1|Block/Mesocosm), data = disp, family = binomial(link = "logit"))
summary(disp_model)
  #Temperature has no significant effect on dispersal