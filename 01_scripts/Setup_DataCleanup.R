library(renv)
library(dplyr)

#Run this to create personal access token with github
#install.packages("gitcreds")
#library(gitcreds)
#gitcreds_set()

#install latex
tinytex::install_tinytex()
tiny::tlmgr_update()


#Reading in an previewing data
df <- read.csv("./00_rawdata/backswimmer_development_dispersal.csv", header =TRUE) 
head(df)
df$Treatment <- as.factor(df$Treatment) #change temperature treatment from numeric to factor
summary(df)

#Creating new column MoltRate that considers how frequently an individual is molting
#Creating new column logMoltRate which log transforms molt rate
df <- df %>%
  mutate(MoltRate = case_when(!is.na(AdulthoodDay) ~ MoltCount/AdulthoodDay, Survival == 0 | is.na(AdulthoodDay)  ~ MoltCount/Death_Day)) %>%
  mutate(logMoltRate = case_when(MoltRate > 0 ~ log(MoltRate)))

  #Rationale for new columns 
  hist(df$MoltRate) #Molt rate is right skewed
  hist(df$logMoltRate) #log transformation helps data fit normal distribution
  
#Exporting revised df to outdata folder
write.csv(df, "./02_outdata/revisedSep18-25_backswimmer_development_dispersal.csv")
  