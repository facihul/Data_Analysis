
StartTme <- Sys.time()
library(ggplot2)
library(e1071)
library(caret) # for confusion matrices

################# load data #####
##########################  Unbalanced sample data ########
#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/NonCodingRNA.RData")
#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/PCodingRNA.RData")

############################# balanced sample data ############
#load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/NonCodingRNA.RData")
#load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/PCodingRNA.RData")

####################### miRNA ##########################
load("/Users/mdfacihulazam/Desktop/ref/RCode/miRNA_seq/new_miRNA_Seq.RData")
##################################################
#data <- PCodingRNA
#data <- nonCodingRNA

############# Class lebels #############
#labels <- c(rep('Class1',500),rep('Class2',462)) # unbalanced lavels
#labels <- c(rep('Class1',212),rep('Class2',212)) # balanced lavels 
#labels <- c(rep('Class1',130),rep('Class2',125)) # unbalanced miRNA
labels <- c(rep('Class1',125),rep('Class2',125)) # balanced miRNA

sampleNum <- length(labels)
row.names(data) <- c(1:sampleNum)
randinds <- sample(1:nrow(data),nrow(data),replace = F) # Random Sample Collection
newLabel <- labels[randinds] # Randomly organize lavels according to sample
data <- data[randinds,]
a <- colnames(data)
# ####################
data <- data.frame(data)


################

data1 <- cbind(data,as.factor(newLabel))
colnames(data1) <- a
data <- data1
num <- ncol(data)

newLabel  <- data[,num]  ##  Data Labels created 
indx <- c(1:length(labels))


################ best PCA set selection #########
g <- c(5,10,15,20,25,30)
r <- c()
p <- c()
allScore <- matrix(nrow = 6,ncol = 6)
num <- ncol(data)
dataP <- data[,-num]
x <- sapply(dataP,as.numeric)
pcadata <- prcomp(x, scale.=TRUE)

for (l in 1:length(g) ){
  
  Pdata <- pcadata$x[,1:g[l]] # 5 ,10,15,20,25,30 PCA
  
  data <- Pdata
  ############  80% traning 20% testdata ##########
  
  start <- 1
  end <- 44
  sample_test <- indx[start:end]   # test index
  sample_train <- indx[-c(start:end)]  # training index
  
  train_in <- data[sample_train,-num] # training sample
  train_out <- newLabel[sample_train] # training  sample label
  
  test_in <- data[sample_test,-num] # test sample
  test_out <- newLabel[sample_test] # test sample
  
  model <- e1071::svm(train_in, train_out, kernel="linear") # SVM model generate
  predsvm1 <- predict(model,test_in) # class prediction

  cm<-confusionMatrix( predsvm1,  test_out,dnn = c("Prediction", "Reference"), 
                       prevalence = NULL,  mode = "sens_spec")  #confusion matrix
  Cmatrix <- cm$table  # confusion matrix
  
  # ######  Two Class ########
  acc <- cm$overall[1]*100 # accuracy
  sensitivity_svm <-  cm$byClass[1] #  sensitivity
  specificity_svm <-  cm$byClass[2] #  specificity
  precision_svm  <-  cm$byClass[5]  #  Precision
  recall_svm  <-  cm$byClass[6]     # recall
  F1 <- cm$byClass[7]               #  F1 Score
  
  r[l] <- sensitivity_svm
  p[l] <- 1-specificity_svm
  allScore[l,]  <-  c(acc, sensitivity_svm, specificity_svm, precision_svm , recall_svm,  F1)
  
  
  
}
########################

### allScore #######


all <- data.frame(PREcision = p, Recall=r,PCA=rep(paste0("PCA", g[1:6])))

EndTime <- Sys.time()
ElapsTime <- EndTime -StartTme

print(ElapsTime)

ggplot( all, aes(x=Recall, y=PREcision)) + geom_point(aes(colour = PCA)) + geom_line() +
  labs(title = "Sensitivity vs. 1 - Specificity",
       x = "Sensitivity ",  
       y = "1 - Specificity" ) 
