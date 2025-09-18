---
output:
  html_document: default
  pdf_document: default
---
# Data dictionary for raw data 

Project: Assessing the effects of temperature on backswimmer development and dispersal

Author: Nicole Regimbal


## See below for definitions of the columns included in the dataset

**Block**
- Block describes  when the juvenile backswimmers were introduced to their temperature treatment. Blocks are staggered but temporally overlapping.
- Block number is between 1 and 5, such that backswimmers in block 1 were introduced first and backswimmers in block 5 were introduced last.

**Treatment**
- The temperature treatment that the backswimmer was exposed to. 
- Temperature treatments are: 22°C, 24°C, 26°C, 28°C, and 30°C.

**Mesocosm**
- Describes the specific tank replicate that a backswimmer was contained in. 
- There are 35 total mesocosms.

**ID**
- Unique identifier for each individual backswimmer. 
- 320 backswimmers were used in the experiment.

**MoltXDay**
- Columns 5 to 9 are 'MoltXDay' such that the X is replaced with a number 1 to 5.
- Backswimmers undergo 5 instars to reach adulthood, so each number represents when the molt to the next instar occurs.
- The value in the rows indicate the day of the experiment that the molt occurs on.

**MoltCount**
- Molt count is the total number of molts a backswimmer undergoes. Backswimmers may die before becoming adults, so this column represents the number of molts the individual survives through.

**AdulthoodDay**
- The day of the experiment that the backswimmer molts into an adult

**Survival**
- Indicates whether a backswimmer lives or dies.
- The values are binary. 1 indicates an individual survives and 0 indicates an individual dies.

**Death_Day**
- If an individual dies, this indicates the day of the experiment that death occurs on. 
- If an individaul survives, the value in the row is NA.

**Thorax_Width1**
- Thorax width (mm) of the adult backswimmer immediately after it develops into an adult.

**Mass1**
- Mass (g) of the adult backswimmer immediately after it develops into an adult.

**Ambient_Date**
- When a backswimmer develops into an adult, it is removed from the temperature treatment to ambient temperature until the dispersal assay.
- This column indicates the date that the backswimmer is moved to ambient temperature.
- Backswimmers stay at ambient temperatures for a minimum of 10 days. 
- This is important to match backswimmers that developed into adults at similar times into the same dispersal assay group.

**Thorax_Width2**
- Thorax width (mm) of the adult backswimmer immediately before the dispersal assay.

**Mass2**
- Mass (g) of the adult backswimmer immediately before the dispersal assay.

**Disp_Date**
- The date that the dispersal assay is conducted for the backswimmer.

**Disp_Count**
- The number of dispersal attempts a backswimmer makes in the dispersal propensity assay.

**Dispersed**
- Whether an individual dispersed in the dispersal ability assay.
- Dispersal is absolute, 0 indicates an individual stayed and 1 indicates an individual dispersed.

