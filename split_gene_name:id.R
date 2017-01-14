


#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/NonZero_RNA_Seq.RData") 
load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/NonZero_RNA_Seq.RData")

genedata <- dataN
genedata_ID <- row.names(dataN) # gene_id  
num <- 1:length(genedata_ID) 
Gene_Id <- vector()
Gene_Name <- vector()
for(i in num)
{
  temp <- unlist(strsplit(genedata_ID[i], "|", fixed = TRUE))
  Gene_Name[i] <- temp[1] 
  Gene_Id[i] <- temp[2] 
}

row.names(genedata) <- Gene_Id 
genedata <- data.frame(genedata)
GeneData <- cbind(Gene_Name,genedata)
GeneData <- data.frame(GeneData)
############### replce sample id with sequnce of number #####
a <- colnames(GeneData)
a2 <- colnames(genedata)
a[2:425] <- c(1:424)  # for newRNA_Seq_data comment this line
a2[1:424] <- c(1:424) #  for newRNA_Seq_data comment this line
#a[2:963] <- c(1:962)  # for balance_newRNA_Seq_data comment this line
#a2[1:962] <- c(1:962) #  for balance_newRNA_Seq_data comment this line
colnames(GeneData) <- a
colnames(genedata) <- a2

#############
path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/"
#path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/"
################ dataset with gene_ID ####
fileout = paste0(path2,"filtered_GID_data.RData")
save(genedata, file =  fileout)
################## dataset with GeneName and ID ######
fileout = paste0(path2,"filtered_GName_ID_data.RData")
save(GeneData, file =  fileout)
