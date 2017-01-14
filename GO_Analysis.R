
library(org.Hs.eg.db)
#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/filtered_GName_ID_data.RData")
load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/filtered_GName_ID_data.RData")
path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/"
#path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/"


########## pathway ######## 
pathways <- as.list(org.Hs.egGO2EG) 
genedata <- data
genelistTotal <- genedata$gene_id

tmp <- unlist(pathways)
nm2 <- intersect(genelistTotal,tmp)
NonCodingGene <- setdiff(genelistTotal,nm2)
ProteinCodingGene <-  nm2


############ Collect samples belongs to Protein coding and Non coding genes ####

PCdata <- genedata[ProteinCodingGene,] 
NCdata <- genedata[NonCodingGene,] 


################ For Non coding data ######
name <- NCdata[,1]
NCdata <- NCdata[,-1]
NCdata <- t(NCdata)
colnames(NCdata) <- name
nonCodingRNA <- data.frame(NCdata)

fileout = paste0(path2,"NonCodingRNA.RData")
save(nonCodingRNA,file = fileout)
############## for protein coding data ######
name <- PCdata[,1]
PCdata <- PCdata[,-1]
PCdata <- t(PCdata)
colnames(PCdata) <- name
PCodingRNA <- data.frame(PCdata)

fileout = paste0(path2,"PCodingRNA.RData")
save(PCodingRNA,file = fileout)


