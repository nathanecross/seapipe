##############################################################################################
################################### FOOOF DATASET ############################################
#Author: Aurore A. Perrault
#Date: Winter 2024
#This script combine fooof files into a single dataset for seapipe
##############################################################################################
rm(list=ls(all=TRUE))
#Package
library(readr) #read data
library(dplyr) #organise df
library(tidyr) #clean df
library(tidyverse)
library(ggplot2) #Visualization

### put all the fooofbands.csv in  the same folder (tips: go in OUT/fooof and CtlF fooofb )
setwd("/Users/m.ela/Desktop/DATA_temp/MISPCBTI/OUT/fooof_csv/") #whare are the fooof csv
#x <- read_csv("P31_B_fooofbands.csv")
##Load files into one vector
files <- list.files(pattern = "*.csv") #look for files starting with same name and ending with.csv and put them in a vector
df_list <- lapply(files,read_csv, skip=0) #read each csv of the vector and convert them into multiples df
size <- nrow(bind_rows(df_list)) #numbers of rows for all subject in total
print (size)

#combine all individual csv into 1 df
require(magicfor)
magic_for() #function to store loop
for(j in 1:length(files)) {
  a <- data.frame(substr(files[j],1,19)) #extract name subject - number depend on your ID coding
  nrows = nrow(data.frame(df_list[j])) #how many segment per subject
  id <- t(data.frame(rep(a,nrows))) #put name subject for each rows
  put(id) #works with magic_for - store variable
}
n <- magic_result_as_vector() #stored variables 
ids <- data.frame(n) #put in a df
magic_free() #stop function magic

data_allsuj <- bind_rows(df_list) #bind all subject .csv
df <- cbind(ids,data_allsuj) #add ID column

## Clean
rm(data_allsuj, df_list, ids, n, size, files) #remove some useless df
clean <- df
clean <- tidyr::separate(clean, n, c("sub", "ses"), sep = "_") 
clean <- select(clean, -...1, -`Start time`, -`End time`,-Duration,  
                -Stitches, -Stage, -Cycle, -`Event type`, -`9-16 Hz`)
clean <- select(clean, -grep("_BW", names(clean)), 
                -grep("_Offset", names(clean)), -grep("_Exponent", names(clean)))

## Biggest PW = peak
peak <- clean 
peak$max <- pmax(peak$`9-16 Hz_peak1_PW`, peak$`9-16 Hz_peak2_PW`,
               peak$`9-16 Hz_peak3_PW`, na.rm = TRUE)

indices <- data.frame(which(select(peak,-("max")) == peak$max, arr.ind = TRUE))
indices$col <- indices$col  - 1
indices <- arrange(indices, row)
mask <- matrix(data=TRUE,nrow=nrow(indices),ncol=1)
for (row in 2:nrow(indices)) {
  if (indices$row[row] == indices$row[row-1]) {
    mask[row] <- FALSE
  }
}

col <- indices$col[mask]

peak$highCF <- NA
for(row in 1:nrow(peak)) {
  peak$highCF[row] <- peak[row,col[row]]
}

## extract .csv and save files
long <- select(peak, sub, ses, Channel, highCF)

## spread long to wide
long$Channel <- ifelse(long$Channel=="C3", "Cz", 
                        ifelse(long$Channel=="Cz", "Cz", 
                               ifelse(long$Channel=="F3", "Fz", 
                                      ifelse(long$Channel=="Fz", "Fz", 
                                             ifelse(long$Channel=="P3", "_REF",
                                                    ifelse(long$Channel=="_REF", "_REF",'x')))))) #


wide <- long %>% 
  spread(Channel, highCF)

unite <- wide %>% 
  unite(col = "CF_chan", c("Fz", "Cz", "_REF"), sep = ",")

final <- wide %>%
  full_join(unite)


setwd("/Users/m.ela/Desktop/DATA_temp/MISPCBTI/OUT/")
write_csv(final, 'sigma_peak_9-16Hz_NREM2-NREM3.csv')


