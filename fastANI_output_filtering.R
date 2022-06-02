##############################################################################################################
#Packages
##############################################################################################################

list.of.packages <- c( "seqinr", "gplots", "reshape2", "stringr", "readr", "dplyr", "BiocManager")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,
                                          repos = "https://cran.uni-muenster.de/")
if(!require(ComplexHeatmap, quietly=TRUE)){
  BiocManager::install("ComplexHeatmap")
}

library(seqinr)
library(reshape2)
library(ComplexHeatmap)
library(gplots)
library(dplyr)
library(stringr)
library(readr)
library(svglite)

##############################################################################################################
#System variables
##############################################################################################################
oDir = Sys.getenv("oDir")
iter = Sys.getenv("iter")
id = Sys.getenv("id")

print(paste("############################"))
print(paste("R Script input parameters:"))
print(paste("output dir =", oDir))
print(paste("iteration =", iter))
print(paste("average nucleotide identity =", id))
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
#import FastANI report
##############################################################################################################
file = file.path(oDir, "05_fastANI", iter, paste0(iter,".fastANI_out.txt"))
fastANI = read_delim(file,
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE)
names(fastANI) = c("query_genome", "reference_genome", "ANI_value", "n_bidir_fragment_mappings", "total_query_mappings")
rm(file)

#remove directory path from sequences so that its easier to work with them
path = file.path(oDir, "04_checkv", iter, "phages/")
path = "/bioinf/home/benedikt.heyerhoff/planktotrons_1/04_checkv/BENG-vm-1_S15_L001/phages/"
#
print(paste("CLEANING FastANI OUTPUT FILE"))
for(i in 1:nrow(fastANI)){
  fastANI[i, "query_genome"] = gsub(paste0(path), "", fastANI[i, "query_genome"])
  fastANI[i, "reference_genome"] = gsub(paste0(path), "", fastANI[i, "reference_genome"])
}
rm(i)


##############################################################################################################
#Heatmap1 of FastANI output
##############################################################################################################
fastANI.matrix = acast(fastANI, query_genome~reference_genome, value.var = "ANI_value")
fastANI.matrix[is.na(fastANI.matrix)] = id

#this whole numeric conversion has to be done because R on the cluster somehow converts the matrix to character...
fastANI.matrix = as.data.frame(fastANI.matrix)
fastANI.matrix[] <- lapply(fastANI.matrix, as.numeric)
####

#create breaks and gradient
breaks = seq(min(fastANI.matrix), max(100), length.out=100)
gradient1 = colorpanel( sum( breaks[-1] <= 95 ), "red", "white" )
gradient2 = colorpanel( sum( breaks[-1] > 95 & breaks[-1] <= 100), "white", "blue" )
#color vector for heatmap
cols = c(gradient1, gradient2)
#path where heatmap will be exported to
path = file.path(oDir, "05_fastANI", iter, "heatmaps")
dir.create(path)

#export as PDF
pdf(paste0(path, "/", "all_contigs_FastANI_Heatmap.pdf"))
  heatmap.2(as.matrix(fastANI.matrix), scale = "none", trace = "none", col = cols, cexRow=.30, cexCol=.30)
dev.off()

rm(fastANI.matrix, cols, gradient1, gradient2, path, breaks)

#################################################################################################################
#filter any entries with identity lower than $id from id input flag in shell script
#################################################################################################################
print(paste0("FILTERING BY ", id, "% AVERAGE NUCLEOTIDE IDENTITY"))
fastANI %>%
  filter(fastANI[, 3] >= id) -> fastANI.id

fastANI.id$query_length = NA
fastANI.id$reference_length = NA
for(i in 1:nrow(fastANI.id)) {
  fastANI.id[i,"query_length"] = as.numeric(gsub(".*;([0-9]+)bp.*", "\\1", fastANI.id[i,"query_genome"]))
  fastANI.id[i,"reference_length"] = as.numeric(gsub(".*;([0-9]+)bp.*", "\\1", fastANI.id[i,"reference_genome"]))
  }
rm(i)

#################################################################################################################
#Import phage sequences
#################################################################################################################
print(paste("IMPORTING CHECKV PHAGE CONTIGS"))
directory = file.path(oDir, "04_checkv", iter, "phages")
contig.list = list.files(path = directory, pattern = ".fna", all.files = TRUE)

phage.sequences = list()
for(i in 1:length(contig.list)) {
  file = file.path(oDir, "04_checkv", iter, "phages", contig.list[i]) 
  tmp = read.fasta(file = file, as.string = TRUE, whole.header = TRUE, 
                               seqtype = "DNA") 
  phage.sequences = append(phage.sequences, tmp)
}
rm(i, file, tmp)

#convert viral sequences to upper case
for(i in 1:length(phage.sequences)){
  phage.sequences[i] = toupper(phage.sequences[i])
}
rm(i)


#################################################################################################################
#quality score expansion
#################################################################################################################
pattern = c("Low-quality", "Medium-quality", "High-quality", "Complete")
  #query
fastANI.id$query_quality=NA
fastANI.id$reference_quality=NA

a = c("query_genome", "reference_genome")
b = c("query_quality", "reference_quality")

for(i in 1:nrow(fastANI.id)){
  for(s in 1:length(a)){
    for(l in 1:length(b)){
      for(k in 1:length(pattern)){
        if (str_detect(fastANI.id[i, a[s]], pattern[k])) {
          fastANI.id[i, b[l]] =  pattern[k]
        }
      }
    }
  }
}
rm(i,k,l,s,a,b)

#add a numerical score to quality in order to make sorting easier
fastANI.id %>%
  mutate(query_q = case_when(
    query_quality == "Low-quality"  ~ 1,
    query_quality == "Medium-quality"  ~ 2,
    query_quality == "High-quality"  ~ 2,
    query_quality == "Complete" ~ 4)) -> fastANI.id

fastANI.id %>%
  mutate(reference_q = case_when(
    reference_quality == "Low-quality"  ~ 1,
    reference_quality == "Medium-quality"  ~ 2,
    reference_quality == "High-quality"  ~ 2,
    reference_quality == "Complete" ~ 4)) -> fastANI.id

#identify duplicated contig entries and keep contigs of higher quality and length
remove.contigs = list()
for(i in 1:nrow(fastANI.id)){
  if (fastANI.id[i, "query_genome"] == fastANI.id[i, "reference_genome"]) {
    next
  }
  if (fastANI.id[i, "query_genome"] != fastANI.id[i, "reference_genome"] & 
      fastANI.id[i, "reference_length"] !=  fastANI.id[i, "query_length"]) {
    if (fastANI.id[i, "query_q"] < fastANI.id[i, "reference_q"] | 
        fastANI.id[i, "query_q"] == fastANI.id[i, "reference_q"] &
        fastANI.id[i, "query_length"] < fastANI.id[i, "reference_length"]) {
      remove.contigs = append(remove.contigs, fastANI.id[i, "query_genome"])
      } else {
      remove.contigs = append(remove.contigs, fastANI.id[i, "reference_genome"])
    } 
  }
}
rm(i)

#unlist to vector
remove.contigs = unlist(remove.contigs, use.names = FALSE)
remove.contigs = remove.contigs[!duplicated(remove.contigs)]
names(phage.sequences) = gsub("(^.*)", "\\1.fna", names(phage.sequences))

#remove duplicate sequences
isNameInIndex <- names(phage.sequences) %in% remove.contigs
duplicated.phage.sequences = phage.sequences[isNameInIndex]
phage.sequences = phage.sequences[!isNameInIndex]

#################################################################################################################
#export phage sequences
#################################################################################################################
dir = file.path(oDir, "05_fastANI", iter, "phages")
dir.create(dir, showWarnings = TRUE)

print(paste("EXPORTING FastANI SORTED PHAGES"))
for(i in 1:length(phage.sequences)) {
  export = phage.sequences[i]
  file = file.path(oDir, "05_fastANI", iter, "phages", names(export))
  write.fasta(export, names(export), file, open = "w")  
}
rm(i)

#################################################################################################################
#export duplicated phages
#################################################################################################################
dir = file.path(oDir, "05_fastANI", iter, "duplicated_phages")
dir.create(dir, showWarnings = TRUE)

print(paste("EXPORTING DUPLICATED PHAGES"))
for(i in 1:length(duplicated.phage.sequences)) {
  export = duplicated.phage.sequences[i]
  file = file.path(oDir, "05_fastANI", iter, "duplicated_phages", names(export))
  write.fasta(export, names(export), file, open = "w")  
}
rm(i)

#################################################################################################################
#Heatmap2 of filtered phages
#################################################################################################################
fastANI.no_id = setdiff(fastANI, fastANI.id[, c(1:5)])
fastANI.matrix <- acast(fastANI.no_id, query_genome~reference_genome, value.var="ANI_value")
fastANI.matrix[is.na(fastANI.matrix)] <- id

#this whole numeric conversion has to be done because R on the cluster somehow converts the matrix to character...
fastANI.matrix = as.data.frame(fastANI.matrix)
fastANI.matrix[] <- lapply(fastANI.matrix, as.numeric)
####

#create breaks and gradient
breaks = seq(min(fastANI.matrix), max(100), length.out=100)
gradient1 = colorpanel( sum( breaks[-1] <= 95 ), "red", "white" )
gradient2 = colorpanel( sum( breaks[-1] > 95 & breaks[-1] <= 100), "white", "blue" )
#color vector for heatmap
cols = c(gradient1, gradient2)
#path where heatmap will be exported to
path = file.path(oDir, "05_fastANI", iter, "heatmaps")
#export as pdf
pdf(paste0(path, "/", "non_identical_contigs_FastANI_Heatmap.pdf"))
  heatmap.2(as.matrix(fastANI.matrix), scale = "none", trace = "none", col = cols, cexRow=.30, cexCol=.30)
dev.off()
rm(fastANI.matrix, cols, gradient1, gradient2, path, breaks)

#################################################################################################################
#Heatmap3 of duplicated phages
#################################################################################################################
fastANI.matrix <- acast(fastANI.id, query_genome~reference_genome, value.var="ANI_value")
fastANI.matrix[is.na(fastANI.matrix)] <- id

#this whole numeric conversion has to be done because R on the cluster somehow converts the matrix to character...
fastANI.matrix = as.data.frame(fastANI.matrix)
fastANI.matrix[] <- lapply(fastANI.matrix, as.numeric)
####

#create breaks and gradient
breaks = seq(min(fastANI.matrix), max(100), length.out=100)
gradient1 = colorpanel( sum( breaks[-1] <= 95 ), "red", "white" )
gradient2 = colorpanel( sum( breaks[-1] > 95 & breaks[-1] <= 100), "white", "blue" )
#color vector for heatmap
cols = c(gradient1, gradient2)
#path where heatmap will be exported to
path = file.path(oDir, "05_fastANI", iter, "heatmaps")
#export as pdf
pdf(paste0(path, "/", "duplicated_contigs_FastANI_Heatmap.pdf"))
  heatmap.2(as.matrix(fastANI.matrix), scale = "none", trace = "none", col = cols, cexRow=.30, cexCol=.30)
dev.off()
rm(fastANI.matrix, cols, gradient1, gradient2, path, breaks)

#clear everything and exit
rm(list = ls())
print(paste("DONE EXPORTING"))
print(paste("EXITING"))

