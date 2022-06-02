##############################################################################################################
#Packages
##############################################################################################################

list.of.packages <- c( "seqinr", "dplyr", "readr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,
                                          repos = "https://cran.uni-muenster.de/")

library(seqinr)
library(dplyr)
library(readr)

##############################################################################################################
#System variables
##############################################################################################################
oDir = Sys.getenv("oDir")
iter = Sys.getenv("iter")

print(paste("############################"))
print(paste("R Script input parameters:"))
print(paste("output dir =", oDir))
print(paste("iteration =", iter))
print(paste("############################"))

##############################################################################################################
#set working directory
##############################################################################################################
directory = file.path(oDir, "02_Virsorter2", iter) #set working directory
setwd(directory)
print(paste("############################"))
print(paste("Setting Working Directory to:"))
print(paste(directory))
print(paste("############################"))

##############################################################################################################
#import virsorter2 results table
##############################################################################################################
file = file.path(oDir, "02_Virsorter2", iter, paste0(iter,"_virsorter_viral_score.tsv"))
virsorter2.tbl = read_delim(file,
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)


##############################################################################################################
#import virsorter2 phages
##############################################################################################################
files = list.files(path = directory, pattern = ".f", all.files = TRUE)
virsorter2.seq = list()
for(i in 1:length(files)) {
  tmp = read.fasta(file = files[i], as.string = TRUE, whole.header = TRUE, seqtype = "DNA") 
  virsorter2.seq = append(virsorter2.seq, tmp)
}
rm(tmp, i)

#convert viral sequences to upper case and 
for(i in 1:length(virsorter2.seq)){
  virsorter2.seq[i] = toupper(virsorter2.seq[i])
}
rm(i)

#order seq list and table to same order
virsorter2.seq = virsorter2.seq[order(names(virsorter2.seq))]
virsorter2.tbl = virsorter2.tbl[order(virsorter2.tbl$seqname), ]

#rename sequences iteratively
for(i in 1:nrow(virsorter2.tbl)){
  seqname = paste0(gsub("\\|", ";", virsorter2.tbl[i, "seqname"]), ";", virsorter2.tbl[i, "max_score_group"], 
                   ";", virsorter2.tbl[i, "length"], "_bp")
  seqname = gsub(";;", ";", seqname)
  names(virsorter2.seq)[i] = seqname
}
rm(seqname, i)

##############################################################################################################
#Export virsorter2 phages
##############################################################################################################
export.dir = file.path(oDir, "03.1_all_viral_contigs")
dir.create(export.dir)
export.dir = file.path(oDir, "03.1_all_viral_contigs", iter)
dir.create(export.dir)

export.dir.seq = file.path(export.dir, paste0(iter,".phage_contigs_virsorter.fna"))

write.fasta(virsorter2.seq, names(virsorter2.seq), export.dir.seq, open = "w")  

#clear everything and exit
rm(list = ls())
print(paste("DONE EXPORTING"))
print(paste("EXITING"))

