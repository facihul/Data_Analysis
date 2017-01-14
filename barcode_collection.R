# This module downloads RNA_Seq data of selected barcode from clinical info . 
# Here TCGAbilinks package is used to download data from TCGA data portal 

###########  clinical data #####

cli_pat <- read.delim(file.choose(), sep = '\t') ## Clinicle data

########################## select barcode of patients of 4 stages ######

a=1
b=1
c=1
d=1
count <- 1:1099;
barcode1A <- vector();
barcode2A <- vector();
barcode2B <- vector();
barcode3A <- vector();

#### barcode Selection according to stages ################

for (i in count) {
  
  if (cli_pat$ajcc_pathologic_tumor_stage[i] == "Stage I")
  {
    barcode1A[a] <- cli_pat$bcr_patient_barcode[i]
    a = a+1;
  }
  else if (cli_pat$ajcc_pathologic_tumor_stage[i] == "Stage IIA")
  {
    barcode2A[b] <- cli_pat$bcr_patient_barcode[i]
    b = b+1;
  }
  else if (cli_pat$ajcc_pathologic_tumor_stage[i] == "Stage IIB")
  {
    barcode2B[c] <- cli_pat$bcr_patient_barcode[i]
    c = c+1;
  }
  else if(cli_pat$ajcc_pathologic_tumor_stage[i] == "Stage IIIA")
  {
    barcode3A[d] <- cli_pat$bcr_patient_barcode[i]
    d = d+1;
  }
}


bar1_mat <- matrix(cli_pat$bcr_patient_barcode[barcode1A])
bar2A_mat <- matrix(cli_pat$bcr_patient_barcode[barcode2A])
bar2B_mat <- matrix(cli_pat$bcr_patient_barcode[barcode2B])
bar3A_mat <- matrix(cli_pat$bcr_patient_barcode[barcode3A])
################ save barcode in a text file if you want ################
# path <-  ............
#filename <- paste0(path,'filename')
#save(bar3A_mat,filename ) 


