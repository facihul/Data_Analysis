
StartTme <- Sys.time()
library(e1071)
library(caret) # for confusion matrices
fold <-4 # 5 fold cross validation
######### date 08.01.2017 modified data ######

################# load data #####
##########################  Unbalanced sample data ########
#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/NonCodingRNA.RData")
#load("/Users/mdfacihulazam/Desktop/ref/RCode/newRNA_Seq_data/PCodingRNA.RData")

############################# balanced sample data ############
#load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/NonCodingRNA.RData")
load("/Users/mdfacihulazam/Desktop/ref/RCode/Balanced_newRNA_Seq_data/PCodingRNA.RData")

####################### miRNA ##########################


#load("/Users/mdfacihulazam/Desktop/ref/RCode/miRNA_seq/new_miRNA_Seq.RData")


##########################

data <- PCodingRNA
#data <- nonCodingRNA

############# Class lebels  for SEER system #############
#labels <- c(rep('Class1',500),rep('Class2',462)) # unbalanced lavels
labels <- c(rep('Class1',212),rep('Class2',212)) # balanced lavels 
#labels <- c(rep('Class1',130),rep('Class2',125)) # unbalanced miRNA
#labels <- c(rep('Class1',125),rep('Class2',125)) # balanced miRNA

sampleNum <- length(labels)
row.names(data) <- c(1:sampleNum)
randinds <- sample(1:nrow(data),nrow(data),replace = F) # Random Sample Collection
newLabel <- labels[randinds] # Randomly organize lavels according to sample
data <- data[randinds,]

############ Apply Variance #########
vr <- apply(data,2, var)
hist(log(vr))
#  # # # ## 10 to 18 for bPCodingRNA
#  # # # ## 14 to 16 for bnonCodingRNA
#  # # # ## 20 to 25 for miRNA
# You can try different combination of log(vr)

k1 <- which(log(vr)<=16)
k2 <- which(log(vr)>=14)
k <- intersect(k1,k2)  # feature Sample collection 

#  ##################### feature variance #####
#  # vars <-apply(data,MARGIN = 2,FUN = var)
#  # inputs <- data[,order(vars,decreasing = T)[1:Ngenes]]
#  #  data1 <- inputs
#  ##########################
data1 <- data[,k]
data <- data1
a <- colnames(data)
data <- data.frame(data)
################

data1 <- cbind(data,as.factor(newLabel))
colnames(data1) <- a
data <- data1
num <- ncol(data)
newLabel  <- data[,num]  ##  Data Labels created 

################ PCA#########
# dataP <- data[,-num]
# xa <- sapply(dataP,as.numeric)
# pcadata <- prcomp(xa, scale.=TRUE)
# Pdata <- pcadata$x[,1:200]
# 
# data <- Pdata
########  n fold cross validation (LOOCV) ########

ste <- floor(sampleNum/fold)  # each fold size of data
indx <- c(1:length(labels))
#indx <- row.names(data1)
j<-1
F1 <- c()
acc <- c()
sensitivity_svm <- c()
specificity_svm <- c()
precision_svm <-  c()
recall_svm <- 
  sample_test <- c()
sample_train <- c()
count <- 1

for(i in 1:fold){
  start <- count
  count <- i*ste
  sample_test <- indx[start:count]   # test index
  sample_train <- indx[-c(start:count)]  # training index
  
  train_in <- data[sample_train,-num] # training sample
  train_out <- newLabel[sample_train] # training  sample label
  
  test_in <- data[sample_test,-num] # test sample
  test_out <- newLabel[sample_test] # test sample
  
  model <- e1071::svm(train_in, train_out, kernel="linear") # SVM model generate
  predsvm1 <- predict(model,test_in) # class prediction
  #tab <- table(predsvm1,test_out) # confusion matrix
  #print(tab)
  cm<-confusionMatrix( predsvm1,  test_out,dnn = c("Prediction", "Reference"), 
                       prevalence = NULL,  mode = "sens_spec")  #confusion matrix
  Cmatrix <- cm$table  # confusion matrix
  acc[j] <- cm$overall[1] # accuracy
  ######  macro avaraging is done here ########
  sensitivity_svm[j] <-  cm$byClass[1] # per class avarage sensitivity
  specificity_svm[j] <-  cm$byClass[2] # per class avarage specificity
  precision_svm[j]  <-  cm$byClass[5]  # per class avarage Precision 
  recall_svm[j]  <-  cm$byClass[6]     #per class avarage recall
  
  F1[j] <- cm$byClass[7]               # per class avarage F1 Score
  
  j<-j+1
}

## Mean  calculation #######
acc_mean <- mean(acc)*100 # mean of accuracy
AF1    <-  mean(F1)       # Avarage F1 score
sensitivity_mean <- mean(sensitivity_svm)*100 # sensitivity
specificity_mean <-  mean(specificity_svm)*100 # specificity
precision_mean <-  mean(precision_svm)*100 # precision
recall_mean  <-  mean(recall_svm)*100  # recall 

## Standard Error calculation ###### 
Sd_acc = sd(acc)/sqrt(length(acc)) # mean Standered error
Sd_sens = sd(sensitivity_svm)/sqrt(length(sensitivity_svm)) # mean Standered error
Sd_spec = sd(specificity_svm)/sqrt(length(specificity_svm)) # mean Standered error

### allScore #######

Afone <-  c(precision_mean ,    recall_mean,  AF1)
allScore <- c(acc_mean, Sd_acc, sensitivity_mean, Sd_sens, specificity_mean, Sd_spec)
## Print Results ######
print("Accuracy,    Sd_acc,  sensitivity,    Sd_sens,   specificity,   Sd_spec  ")
print(allScore)
# print("precision_mean, recall_mean,  AF1")  
# print(Afone)

EndTime <- Sys.time()
ElapsTime <- EndTime -StartTme

print(ElapsTime)
