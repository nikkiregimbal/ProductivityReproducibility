library(renv)

#Run this to create personal access token with github
#install.packages("gitcreds")
#library(gitcreds)
#gitcreds_set()


#install latex
tinytex::install_tinytex()
tiny::tlmgr_update()

# Acquiring data and saving to .csv ---------------------------------------
install.packages("lterdatasampler")
library(lterdatasampler) #this contains dataset options

df <- and_vertebrates #choosing this dataset and making it an object
write.csv(df, ".00_rawdata/and_vertebrates.csv") #save dataset in raw data folder
#now the data is locally stored

# Initial data exploration ------------------------------------------------

df <- read.csv("./00_rawdata/and_vertebrates.csv")
head(df)



