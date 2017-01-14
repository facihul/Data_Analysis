
## This module reads only normalized gene data from Downloaded IlluminaHiSeq RNA data ####


path <- "/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/miRNA/data_non_3ABRCA/genedata/" # source path
path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/miRNA/data_non_3ABRCA/newdata/"# destination path  

#file_list = list.files(path = path, pattern="*.rsem.genes.results") # list all the files RNA_seq
file_list = list.files(path = path, pattern="*.mirna.quantification.txt") # list all the files for miRNA

st <- vector()
for (j in 1:length(file_list)){
  #st[j] <- strsplit(file_list[j],".rsem.genes.results") # separate filename as  patient id 
  st[j] <- strsplit(file_list[j],".mirna.quantification.txt")
  
}

Pre_filtered_data <- matrix()  # empty matrix for all patients data 
filein = paste0(path,file_list[1])
temp_data <- read.delim(filein, sep = '\t') # read single file in a matrix. 
Pre_filtered_data <-temp_data
#Pre_filtered_data <- Pre_filtered_data[,-c(3,4)] # take only raw-count RNA_seq
Pre_filtered_data <- Pre_filtered_data[,-c(3,4)] # take only read-count miRNA
for (i in 2:length(file_list)){
  
  filein = paste0(path,file_list[i])
  temp_data <- read.delim(filein, sep = '\t')
  a <- c(temp_data$read_count) # miRNA 
  #a <- c(temp_data$raw_count)
  Pre_filtered_data <- cbind(Pre_filtered_data,a)
  
}

## Rows of the matrix define geneid and column represents sample name / patient id
for (i in 2:ncol(Pre_filtered_data)) {
  colnames(Pre_filtered_data)[i] <- st[i-1] # replace colnames with patient id
}
a <- unique(colnames(Pre_filtered_data))  # check any duplicate exist
newdata <- Pre_filtered_data[,a] 
# 
# #  save the file ############
 fileout = paste0(path2,"newdata.RData")
 save(newdata,file = fileout)

