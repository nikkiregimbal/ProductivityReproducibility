library(ggplot2)
library(dplyr)
library(vegan)
library(lme4)
library(partR2)

df <- read.csv("./02_outdata/revisedSep18-25_backswimmer_development_dispersal.csv")


# Temperature on Dispersal ------------------------------------------------
#Exploring the data
disp <-subset(df, Disp_Count >= 0) #only 30 individuals survived to the dispersal assay

disp_count <- table(disp$Treatment, disp$Disp_Count) #table of individuals in dispersal propensity assay by temperature treatment

sum(disp$Disp_Count) #only 3 individuals attempted dispersal in the dispersal propensity assay

disp_final <- table(disp$Treatment, disp$Dispersed) #table of individuals that dispersed in dispersal ability assay by temperature

sum(disp$Dispersed) #only 1 individual actually dispersed away from the experiment

#Takeaway - dispersal was extremely low in this experiment and the sample size is very small
  #Any analysis lacks power
 
#Since the response is count data, I am conducting a GLMM with poisson distribution
dispcount_model <- glmer(Disp_Count ~ Treatment + (1|Block/Mesocosm), data = disp, family = poisson)
summary(dispcount_model)
  #Temperature has no significant effect on number of dispersal attempts

#Absolute dispersal is binary so I am conducting GLMM with binomial distribution 
disp_model <- glmer(Dispersed ~ Treatment + (1|Block/Mesocosm), data = disp, family = binomial(link = "logit"))
summary(disp_model)
  #Temperature has no significant effect on dispersal