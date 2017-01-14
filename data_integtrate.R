
## data loading  RNA_Seq   ##
# load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/data1ABRCA/newdata/newdata.RData")
# data1 <- newdata # stage 1 - 106 samples
# load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/data2ABRCA/newdata/newdata.RData")
# data2 <- newdata # stage 2A - 394 samples
# load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/data2BBRCA/newdata/newdata.RData")
# data3 <- newdata # stage 2B - 291 samples
# load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/data3ABRCA/newdata/newdata.RData")
# data4 <- newdata # stage 3A -172 samples
######################### loading miRNA data #####################
load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/miRNA/data_non_1BRCA/newdata/newdata.RData")
data1 <- newdata # stage 1 - 15 samples
load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/miRNA/data_non_2ABRCA/newdata/newdata.RData")
data2 <- newdata # stage 2A - 117 samples
load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/miRNA/data_non_2BBRCA/newdata/newdata.RData")
data3 <- newdata # stage 2B - 78 samples
load("/Users/mdfacihulazam/Desktop/ref/RCode/allSeq/miRNA/data_non_3ABRCA/newdata/newdata.RData")
data4 <- newdata # stage 3A -49 samples
############################ check prograssive cancer patient #####
a1 <-colnames(data1) 
a2 <-colnames(data2) 
a3 <-colnames(data3) 
a4 <-colnames(data4) 


###############  equal number of sample data ########
data1A <- data1
data2A <- data2[,c(2:length(a))]
data2B <- data3[,c(2:length(a))]
data3A <- data4[,c(2:length(a))]

 dataA <- cbind(data1A,data2A,data2B,data3A)
 dataT <-  dataA
 row.names(dataT) <- dataT[,1] # assign rowname as gene_id
 data <- dataT[,-1]
 path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/"
############################

############ combine all stages of data in one file ############## 
####### Total 424 samples ##############
# data <- cbind(data1,data2,data3,data4) 
# genename_Dup <- which(colnames(data)=="miRNA_ID") # find duplicate columnname as gene_id/miRNA_ID
# data <- data[,-genename_Dup[2:4]] #  remove duplicate columnname as gene_id
# dataT <-  data
# row.names(dataT) <- dataT[,1] # assign rowname as gene_id
# data <- dataT[,-1]
################# check unique sample ######
a <- unique(colnames(data))
#data <- data[,a]
################# save the data #########

#path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/" # allRNA_Seq data 
path2 <- "/Users/mdfacihulazam/Desktop/ref/RCode/miRNA_Seq/" #  miRNA data
fileout = paste0(path2,"miRNA_Seq.RData")
save(data, file =  fileout)

################## End ###############


