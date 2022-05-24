if(!require('dplyr'))install.packages('dplyr')
if(!require('seqinr'))install.packages('seqinr')
if(!require('stringr'))install.packages('stringr')

library(dplyr)
library(seqinr)
library(stringr)


oDir = Sys.getenv("oDir")
iter = Sys.getenv("iter")
minlength = Sys.getenv("minlength") #parameter for viral contig length

print(paste("############################"))
print(paste("R Script input parameters:"))
print(paste("output dir =", oDir))
print(paste("iteration =", iter))
print(paste("minlength =", minlength))
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
directory = file.path(oDir, "04_checkv", iter, "quality_summary.tsv")#make directory path 
quality = read.delim(directory, fill = TRUE)
rm(directory)#remove directory variable

#filter checkv quality summary by various parameters
quality %>% 
  dplyr::filter(checkv_quality == "Complete") -> complete_phages

quality %>% 
  dplyr::filter(checkv_quality == "High-quality" & contig_length >= minlength) -> high_quality

quality %>% 
  dplyr::filter(checkv_quality == "Medium-quality" & contig_length >= minlength) -> medium_quality

quality %>% 
  dplyr::filter(checkv_quality == "Low-quality" & contig_length >= minlength) -> low_quality

quality %>% 
  dplyr::filter(checkv_quality == "High-quality" & contig_length >= minlength) -> high_quality

quality %>% 
  dplyr::filter(provirus == "Yes" & contig_length >= minlength) -> proviruses


#make a vector contianing all names of dfs
dfs = c("complete_phages", "high_quality", "medium_quality", "low_quality", 
        "complete_phages", "proviruses")
#Place DFs into list
qc_list = list(complete_phages, high_quality, medium_quality, low_quality, complete_phages, proviruses)
names(qc_list) = dfs

###Import phages sequences###
#make directory path for fasta import
print(paste("IMPORTING FASTA FILE(S) WITH VIRAL SEQUENCES"))
directory = file.path(oDir, "04_checkv", iter, "viruses.fna")
viral_sequences = read.fasta(file = directory, as.string = TRUE, whole.header = TRUE, 
                       seqtype = "DNA") #import phage fasta file
directory = file.path(oDir, "04_checkv", iter, "proviruses.fna")
provirus_sequences = read.fasta(file = directory, as.string = TRUE, whole.header = TRUE, 
                                seqtype = "DNA") #import prophage fasta file
rm(directory)#remove directory variable

#convert viral sequences to upper case
for(i in 1:length(viral_sequences)){
  viral_sequences[i] = toupper(viral_sequences[i])
}
rm(i)
#convert proviral sequences to upper case
for(i in 1:length(provirus_sequences)){
  provirus_sequences[i] = toupper(provirus_sequences[i])
}
rm(i)


print(paste("BINNING SEQUENCES"))
#select list items in df "sequences" conditionally based on above determined filtering conditions
#Place DFs into list
complete_phages_seq = viral_sequences[names(viral_sequences) %in% complete_phages$contig_id] 
high_quality_seq = viral_sequences[names(viral_sequences) %in% high_quality$contig_id] 
medium_quality_seq = viral_sequences[names(viral_sequences) %in% medium_quality$contig_id] 
low_quality_seq = viral_sequences[names(viral_sequences) %in% low_quality$contig_id]
provirus_seq = provirus_sequences[names(provirus_sequences)]

#make list of lists
qc_list_seq = list(complete_phages_seq, high_quality_seq, medium_quality_seq, low_quality_seq, provirus_seq)
seq_dfs = c("complete_phages", "high_quality", "medium_quality", "low_quality", "proviruses")
names(qc_list_seq) = seq_dfs


#export these dataframes as individial tables
print(paste("EXPORTING QUALITY TABLES"))
directory = file.path(oDir,"04_checkv", iter, "qc_tables")#make export path
dir.create(directory, showWarnings = TRUE)
for (i in 1:length(dfs)){
  file = file.path(directory, paste0(iter, "_", dfs[i],".tsv"))
  if(nrow(qc_list[[i]]) >= 1) {#only export if there is something in the list
    write.table(qc_list[i], file, sep = ";", row.names = FALSE)
  }
}
rm(directory, i)#remove directory variable


print(paste("EXPORTING PER GENOME FASTA FILES"))
#creating export directory
directory = file.path(oDir,"04_checkv", iter, "phages")#make export path
dir.create(directory, showWarnings = TRUE)

for (i in 1:length(seq_dfs)){
  if (lengths(qc_list_seq[i]) >= 1) {#only export if there is something in the list
    names(qc_list_seq[[i]]) = paste(names(qc_list_seq[[i]]), "checkv", seq_dfs[i], sep = "_")
    for(k in 1:length(qc_list_seq[[i]])){
      export = qc_list_seq[[i]][k]
      names(export) = gsub("/", "|", names(export))
      file = file.path(directory, print(paste0(names(export), ".fna")))
      write.fasta(export, names(export), file, open = "w")
    }
  }
}
#clear everything
rm(list = ls())

print(paste("DONE EXPORTING"))
print(paste("EXITING"))

