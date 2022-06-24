##############################################################################################################
#Packages
##############################################################################################################
list.of.packages <- c( "seqinr", "dplyr", "stringr", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,
                                          repos = "https://cran.uni-muenster.de/")

library(seqinr)
library(dplyr)
library(stringr)
library(tidyr)

##############################################################################################################
#System variables
##############################################################################################################
oDir = Sys.getenv("oDir")
iter = Sys.getenv("iter")
minlength = Sys.getenv("minlength") #parameter for viral contig length

print(paste("############################"))
print(paste("R Script input parameters:"))
print(paste("output dir =", oDir))
print(paste("iteration =", iter))
print(paste("minlength =", minlength))
print(paste("############################"))

##############################################################################################################
#set working directory
##############################################################################################################
directory = file.path(oDir, "04_checkv", iter) #set working directory
setwd(directory)
print(paste("############################"))
print(paste("Setting Working Directory to:"))
print(paste(directory))
print(paste("############################"))
rm(directory)

##############################################################################################################
#import checkv quality report
##############################################################################################################
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
  dplyr::filter(provirus == "Yes" & contig_length >= minlength) -> proviruses


#make a vector contianing all names of dfs
dfs = c("complete_phages", "high_quality", "medium_quality", "low_quality", "proviruses")
#Place DFs into list
qc_list = list(complete_phages, high_quality, medium_quality, low_quality, proviruses)
names(qc_list) = dfs


##############################################################################################################
#Import phages sequences
##############################################################################################################
#Phages
print(paste("IMPORTING FASTA FILE(S) WITH VIRAL SEQUENCES"))
directory = file.path(oDir, "04_checkv", iter, "viruses.fna")
if (file.size(directory) > 0){
  viral_sequences = read.fasta(file = directory, as.string = TRUE, whole.header = TRUE, 
                             seqtype = "DNA") #import phage fasta file
  
  #convert viral sequences to upper case
  for(i in 1:length(viral_sequences)){
    viral_sequences[i] = toupper(viral_sequences[i])
  }
  rm(i)
}else{
  print(paste("No phages in this sample"))
}
#Prophages
directory = file.path(oDir, "04_checkv", iter, "proviruses.fna")
if (file.size(directory) > 0){
  provirus_sequences = read.fasta(file = directory, as.string = TRUE, whole.header = TRUE, 
                                  seqtype = "DNA") #import prophage fasta file
  
  #convert proviral sequences to upper case
  for(i in 1:length(provirus_sequences)){
    provirus_sequences[i] = toupper(provirus_sequences[i])
  }
  rm(i)
}else{
  print(paste("No prophages in this sample"))
}
rm(directory)#remove directory variable

##############################################################################################################
#BINNING
##############################################################################################################
print(paste("BINNING SEQUENCES"))
#select list items in df "sequences" conditionally based on above determined filtering conditions
#Place DFs into list
complete_phages_seq = viral_sequences[names(viral_sequences) %in% complete_phages$contig_id] 
high_quality_seq = viral_sequences[names(viral_sequences) %in% high_quality$contig_id] 
medium_quality_seq = viral_sequences[names(viral_sequences) %in% medium_quality$contig_id] 
low_quality_seq = viral_sequences[names(viral_sequences) %in% low_quality$contig_id]
if(exists("provirus_seq")){
  provirus_seq = provirus_sequences[names(provirus_sequences)]
}

#make list of lists
if(exists("provirus_seq")){
  qc_list_seq = list(complete_phages_seq, high_quality_seq, medium_quality_seq, low_quality_seq, provirus_seq)
  seq_dfs = c("complete_phages", "high_quality", "medium_quality", "low_quality", "proviruses")
  names(qc_list_seq) = seq_dfs
  rm(complete_phages_seq, high_quality_seq, medium_quality_seq, low_quality_seq, provirus_seq)
}else{
  qc_list_seq = list(complete_phages_seq, high_quality_seq, medium_quality_seq, low_quality_seq)
  seq_dfs = c("complete_phages", "high_quality", "medium_quality", "low_quality")
  names(qc_list_seq) = seq_dfs
  rm(complete_phages_seq, high_quality_seq, medium_quality_seq, low_quality_seq)
}

##############################################################################################################
#EXPORt QUALITY TABLES
##############################################################################################################
#export these dataframes as individial tables
print(paste("EXPORTING QUALITY TABLES"))
directory = file.path(oDir,"04_checkv", iter, "qc_tables")#make export path
dir.create(directory, showWarnings = TRUE)
for (i in 1:length(dfs)){
  file = file.path(directory, paste0(iter, "_", dfs[i],".csv"))
  if(nrow(qc_list[[i]]) >= 1) {#only export if there is something in the list
    write.table(qc_list[i], file, sep = ";", row.names = FALSE)
  }
}
rm(directory, i)#remove directory variable

print(paste0("####################################################"))
print(paste0("If a warning message occurs here, it can probably be ignored"))
#Make a quick summary df
df <-
  tidyr::separate(
    data = quality,
    col = contig_id,
    sep = ";",
    into = c("contig", "tool"),
    remove = TRUE
  )
print(paste0("If a warning message occurs here concerning additional pieces discarded in n rows rows, it can be ignored"))
print(paste0("####################################################"))

for(i in 1:nrow(df)) {
  if (str_detect(df[i,"tool"], "full") == TRUE) {
    df[i,"tool"] = gsub("full.*", "virsorter2", df[i,"tool"])
  }
  if (str_detect(df[i,"tool"], "partial") == TRUE) {
    df[i,"tool"] = gsub("partial.*", "virsorter2", df[i,"tool"])
  }
  if (str_detect(df[i,"tool"], "lt2gene") == TRUE) {
    df[i,"lt2gene"] = "yes"
    df[i,"tool"] = gsub(".*", "virsorter2", df[i,"tool"])
  }
}
rm(i)


##############################################################################################################
#create quick summary
##############################################################################################################
df %>%
  group_by(tool, checkv_quality, provirus, {if("lt2gene" %in% names(.)) lt2gene else NULL}) %>% 
  summarise(contig_count = n(),
            Mean_contig_length = mean(contig_length, na.rm=TRUE), 
            Median_contig_length = median(contig_length, na.rm=TRUE)) -> a 

if ("lt2gene" %in% names(a)) {
  a$lt2gene[is.na(a$lt2gene)] = "no"
}


a$handling_by_filter_script = NA
for(i in 1:nrow(a)){
  if (str_detect(a[i, "checkv_quality"], "Not-determined") == TRUE){
    a[i, "handling_by_filter_script"] = "removed by filtering"
  }else{
    if("lt2gene" %in% names(a)){
      if (str_detect(a[i, "lt2gene"], "yes") == TRUE){
      a[i, "handling_by_filter_script"] = "removed by filtering"
      }
      }else{
      a[i, "handling_by_filter_script"] = "retained"
    }
  }
}
print(paste("EXPORTING QUICK SUMMARY"))
directory = file.path(oDir,"04_checkv", iter, paste0(iter,"_FILTER_QUICKSUMMARY.csv"))
write.table(a, directory, sep = ";", row.names = FALSE)


##############################################################################################################
#EXPORT PER GENOME/FRAGMENT FASTA FILES
##############################################################################################################
print(paste("EXPORTING PER GENOME/FRAGMENT FASTA FILES"))
#creating export directory
directory = file.path(oDir,"04_checkv", iter, "phages")#make export path
dir.create(directory, showWarnings = TRUE)

for (i in 1:length(seq_dfs)){
  if (lengths(qc_list_seq[i]) >= 1) {#only export if there is something in the list
    for(k in 1:length(qc_list_seq[[i]])){
      export = qc_list_seq[[i]][k]
      if (str_detect(names(export), "lt2gene") == TRUE) {
        print(paste("SKIPPING THIS SEQUENCE AS IT CONTAINS LESS THAN 2 HALLMARK GENES"))
        next
      }
      
      if (str_detect(names(export), "full") == TRUE) {
        #names(export) = gsub("(virsorter2).*", "\\1", names(export))
        n = gsub("_[0-9] .*", "", names(export)) #this is in case the sequence is a prophage sequence. it cuts the name 
        n = paste0(gsub("full", "virsorter2",
                        gsub(";", "_",
                                  gsub("([0-9]+)_bp", "", n))),
                   ifelse(quality[which(quality == n), "provirus"] == "Yes", 
                          paste0("_", quality[which(quality == n), "checkv_quality"], "_provirus", "_",
                                 quality[which(quality == n), "contig_length"], "bp"), 
                          paste0("_", quality[which(quality == n), "checkv_quality"],"_", 
                                  quality[which(quality == n), "contig_length"], "bp")))
        #find the name of current iteration in quality dataframe and add contig quality, length and if
        #complete or partial to name so it can be exported as fasta header
        names(export) = n
      }
      
      if (str_detect(names(export), "partial") == TRUE) {
        n = gsub("_[0-9] .*", "", names(export))
        n = paste0(gsub("full", "virsorter2",
                        gsub(";", "_",
                             gsub("([0-9]+)_bp", "", n))),
                   ifelse(quality[which(quality == n), "provirus"] == "Yes", 
                          paste0("_", quality[which(quality == n), "checkv_quality"], "_provirus", "_",
                                 quality[which(quality == n), "contig_length"], "bp"), 
                          paste0("_", quality[which(quality == n), "checkv_quality"],"_", 
                                 quality[which(quality == n), "contig_length"], "bp")))
        #find the name of current iteration in quality dataframe and add contig quality, length and if
        #complete or partial to name so it can be exported as fasta header
        names(export) = n
      }
      
      if (str_detect(names(export), "vibrant") == TRUE) {
        n = gsub("_[0-9] .*", "", names(export))
        n = paste0(gsub(";", "_", n), 
                   ifelse(quality[which(quality == n), "provirus"] == "Yes", 
                          paste0("_", quality[which(quality == n), "checkv_quality"], "_provirus", "_",
                                 quality[which(quality == n), "contig_length"], "bp"), 
                          paste0("_", quality[which(quality == n), "checkv_quality"], "_",
                                 quality[which(quality == n), "contig_length"], "bp")))
        #find the name of current iteration in quality dataframe and add contig quality, length and if
        #complete or partial to name so it can be exported as fasta header
        names(export) = n
      }
      
      if (str_detect(names(export), "metaviralspades") == TRUE) {
        n = gsub("_[0-9] .*", "", names(export))
        n = paste0(gsub(";", "_", n),  
                   ifelse(quality[which(quality == n), "provirus"] == "Yes", 
                          paste0("_", quality[which(quality == n), "checkv_quality"], "_provirus", "_",
                                 quality[which(quality == n), "contig_length"], "bp"), 
                          paste0("_", quality[which(quality == n), "checkv_quality"], "_", 
                                 quality[which(quality == n), "contig_length"], "bp")))
        #find the name of current iteration in quality dataframe and add contig quality, length and if
        #complete or partial to name so it can be exported as fasta header
        names(export) = n
      }
      file = file.path(directory, print(paste0(names(export), ".fna")))
      write.fasta(export, names(export), file, open = "w")
    }
  }
}
#clear everything
rm(list = ls())

print(paste("DONE EXPORTING"))
print(paste("EXITING"))

