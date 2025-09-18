---
output:
  html_document: default
  pdf_document: default
---
# Data dictionary for out data 

Project: Assessing the effects of temperature on backswimmer development and dispersal

Author: Nicole Regimbal

#### **See data dictionary in raw data folder for definitons of all original columns in the dataset. Below defines only new columns.**

**MoltRate**
- Describes the frequency at which backswimmers molted,a metrix for developmental rate.
- Calculated as the number of molts / days it took to reach adulthood.
- If an individual dies before becoming an adult, molt rate is calcualted as number of molts / days lived.

**logMoltRate**
- log10 transformation fo MoltRate column. 
- This changes the distribution of data from right-skewed to gaussian. 