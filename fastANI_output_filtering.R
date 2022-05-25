if(!require('dplyr'))install.packages('dplyr')
if(!require('seqinr'))install.packages('seqinr')
if(!require('stringr'))install.packages('stringr')
if(!require('readr'))install.packages('readr')

library(dplyr)
library(seqinr)
library(stringr)
library(readr)


oDir = Sys.getenv("oDir")
iter = Sys.getenv("iter")

print(paste("############################"))
print(paste("R Script input parameters:"))
print(paste("output dir =", oDir))
print(paste("iteration =", iter))
print(paste("############################"))

#set working directory
directory = file.path(oDir, "04_checkv", iter) #set working directory
setwd(directory)
print(paste("############################"))
print(paste("Setting Working Directory to:"))
print(paste(directory))
print(paste("############################"))
rm(directory)

#import checkv quality report
directory = file.path(oDir, "05_fastANI", iter, paste0(iter,".fastANI_out.txt"))
fastANI = read_delim("~/Desktop/05_fastANI/BENG-vm-1_S15_L001/BENG-vm-1_S15_L001.fastANI_out.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE)
rm(directory)

#remove path from sequencees so that its easier to work with them
path = file.path(oDir, "04_checkv", iter, "phages/")
path = "/bioinf/home/benedikt.heyerhoff/planktotrons_1/04_checkv/BENG-vm-1_S15_L001/phages/"
for(i in 1:nrow(fastANI)){
  fastANI[i, 1] = gsub(paste0(path), "", fastANI[i, 1])
  fastANI[i, 2] = gsub(paste0(path), "", fastANI[i, 2])
}

rm(i)





















#clear everything
rm(list = ls())

print(paste("DONE EXPORTING"))
print(paste("EXITING"))

