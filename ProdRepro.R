library(renv)

#Run this to create personal access token with github
#install.packages("gitcreds")
#library(gitcreds)
#gitcreds_set()

library(lterdatasampler) #this contains dataset options
#write.csv(df, "and_vertebrates.csv") saving chosen dataset

df <- read.csv("and_vertebrates.csv")

#install latex
tinytex::install_tinytex()
tiny::tlmgr_update()
