
load("/Users/mdfacihulazam/Desktop/ref/RCode/miRNA_seq/miRNA_Seq.RData")
#load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/RNA_Seq.RData")
#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/RNA_Seq.RData")
data <- t(data)
n <- ncol(data)
co <- c()
count <- 1
for(i in 1:n){
  temp <- data[,i]
  a <- which(temp > 0) ## to find above 0 value gene
  if (length(a)==0) co[count] <-i
  else co[count] <- 0
  count <- count+1
}
b <- which(co!=0) # find index conteins 0 value
datatest <- data[,-b] # remove zero value 
#dataN <- t(datatest)
data <- data.frame(datatest)
#path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/"
#path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/"
path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/miRNA_seq/" # for miRNA 

#fileout = paste0(path2,"NonZero_RNA_Seq.RData") # for all RNA
fileout = paste0(path2,"new_miRNA_Seq.RData")  # for miRNA  
save(data, file =  fileout)

