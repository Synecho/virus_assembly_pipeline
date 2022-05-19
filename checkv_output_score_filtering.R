library(dplyr)
library(seqinr)

oDir = Sys.getenv("oDir")
currentDir = Sys.getenv("currentDir")
iter = Sys.getenv("iter")
compl = Sys.getenv("compl") #parameter for checkv completeness
minlength = Sys.getenv("minlength") #parameter for viral contig length

print(paste("############################"))
print(paste("R Script input parameters:"))
print(paste("output dir:", oDir))
print(paste("current dir:", currentDir))
print(paste("iteration:", iter))
print(paste("completeness:", compl))
print(paste("minlength:", minlength))
print(paste("############################"))

#import checkv quality report
directory = file.path(currentDir, "quality_summary.tsv")#make directory path 
quality = read.delim(directory, fill = TRUE)
rm(directory)#remove directory variable

#filter checkv quality summary by various parameters
q_filtered = dplyr::filter(quality, completeness > compl)
length_filtered = dplyr::filter(quality, contig_length > minlength)
high_quality =  dplyr::filter(quality, checkv_quality == "High-quality")
med_quality = dplyr::filter(quality, checkv_quality == "Medium-quality")
low_quality = dplyr::filter(quality, checkv_quality == "Low-quality")
high_score = dplyr::filter(quality, checkv_quality == "High-quality" && quality$checkv_quality >= 80)
complete_phages = dplyr::filter(quality, checkv_quality == "Complete")


#make a vector contianing all names of dfs
dfs = c("q_filtered", "length_filtered", "high_quality", "med_quality", "low_quality", 
        "high_score", "complete_phages")
#Place DFs into list
qc_list = list(q_filtered, length_filtered, high_quality, med_quality, low_quality, 
               high_score, complete_phages)
names(qc_list) = dfs

###Import phages sequences###
#make directory path for fasta import
print(paste("IMPORTING FASTA FILE WITH (VIBRANT) VIRAL SEQUENCES"))
directory = file.path(currentDir, "viruses.fna")
sequences = read.fasta(file = directory, as.string = TRUE, whole.header = TRUE, seqtype = "DNA")#import phage fasta file
rm(directory)#remove directory variable

#convert sequences to upper case
for(i in 1:length(sequences)){
  sequences[i] = toupper(sequences[i])
}
rm(i)


print(paste("FILTERING SEQUENCES"))
#select list items in df "sequences" conditionally based on above determined filtering conditions
q_filtered_seq = sequences[q_filtered$contig_id]
length_filtered_seq = sequences[length_filtered$contig_id]
high_quality_seq = sequences[high_quality$contig_id]
med_quality_seq = sequences[med_quality$contig_id]
low_quality_seq = sequences[low_quality$contig_id]
high_score_seq = sequences[high_score$contig_id]
complete_phages_seq = sequences[complete_phages$contig_id]

#make a vector contianing all names of seq_dfs
seq_dfs = c("q_filtered_contigs", "length_filtered_contigs", "high_quality_contigs", "med_quality_contigs", 
            "low_quality_contigs", "high_score_contigs", "complete_phages_contigs")
#Place DFs into list
qc_list_seq = list(q_filtered_seq, length_filtered_seq, high_quality_seq, med_quality_seq, 
                   low_quality_seq, high_score_seq, complete_phages_seq)
names(qc_list_seq) = seq_dfs

rm(q_filtered_seq, length_filtered_seq, high_quality_seq, med_quality_seq, 
   low_quality_seq, high_score_seq, complete_phages_seq)

#export these dataframes as individial tables
directory = file.path(oDir,"checkv", "qc_tables")#make export path

for (i in 1:length(dfs)){
  file = file.path(directory, paste0(iter, "_", dfs[i],".tsv"))
  if(nrow(qc_list[[i]]) >= 1) {#only export if there is something in the list
    write.table(qc_list[i], file, sep = ";", row.names = FALSE)
  }
}
rm(directory, i)#remove directory variable

#export filtered fastas as individual files WHILE ALSO REMOVE NON-VIRAL CONTIGS
print(paste("SAVING FILTERED FASTA FILES"))
print(paste("AND"))
print(paste("REMOVING NON-VIRAL CONTIGS"))
for (i in 1:length(seq_dfs)){
  directory = file.path(oDir,"checkv", seq_dfs[i])#make export path
  file = file.path(directory, paste0(iter, ".fna"))
  if (lengths(qc_list_seq[i]) >= 1) {#only export if there is something in the list
    ls.new = qc_list_seq[[i]]
    #ls.new = ls[-which(non_viral_contigs$contig_id %in% names(ls))]
    write.fasta(ls.new, names(qc_list_seq[[i]]), file, open = "w")
  }
}
rm(directory, file, i, ls.new)#remove directory variable

#Output for VirHostMatcher-Net
ls.new = qc_list_seq$length_filtered
#ls.new = ls[-which(non_viral_contigs$contig_id %in% names(ls))]
#exporting only viral contigs
for(i in 1:length(ls.new)){
  if (lengths(ls.new[i]) >= 1) {#only export if there is something in the list
    directory = file.path(oDir, "checkv", "VirHostMatcher", iter) #make export path
    file = file.path(directory, print(paste0(names(ls.new[i]), ".fna")))
    write.fasta(ls.new[i], names(ls.new[i]), file, nbchar = 60, open = "w")
  }
}
rm(ls.new, i, file, directory)
rm(list = ls(all.names = TRUE)) #clear all objects includes hidden objects






