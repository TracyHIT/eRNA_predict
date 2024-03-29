###########Train and test model #############
#' First train in 5-fold and validated to select the best parameter(mtry,ntree).
#' Then train the model with the best parameter and test the model 5 times
#' output: RF model saved in Result/Model/, error_matrix saved in Result/Summary/
#'
#' @param data_file: the matrix Rdata generated by make_model_feature_index, last colomn is Y
#' @param cell: the cell type name
#' @param Feature_index: the feature index, please refer to make_model_feature_index
RF_train_self_tuning_fuzz<- function(data_file,cell,Feature_index=3)
{
  library("randomForest")
  library("pROC")
  library("caret")
  library("doParallel")
  library("foreach")
  load(data_file)
  source("Script/make_model_feature_index.R")
  accepts <- Freature_index_make(accepts,Feature_index)
  accepts <- accepts[which(accepts$Ac_Y ==0 | accepts$Ac_Y ==1),]
  ######5-fold
  len <- nrow(accepts)
  k = 5
  mysplit = function(k,len){
    pool = c(1:len)
    seg = as.integer(len/k)
    train = as.data.frame(matrix(nrow = (len - seg)))
    test = as.data.frame(matrix(nrow = seg))
    for (i in 1 : k){
      ctest = sample(pool,seg,replace = FALSE)
      train[i] = setdiff(c(1:len),ctest)
      test[i] = ctest
      pool = setdiff(pool,ctest)
    }
    out = list(one = train, two = test)
    return(out)
  }
  t_mysplit<- mysplit(k,len)
  train_index = t_mysplit$one
  test_index = t_mysplit$two
  
  


  cl <- makeCluster(4)
  registerDoParallel(cl)
  result <- foreach(tmp_mtry = rep(c(2,4,6,8),5),
                    tmp_ntree =  c(rep(100,4),rep(200,4),rep(300,4),rep(400,4),rep(500,4))) %dopar% {
                      library("randomForest")
                      OOB_err_temp <- NULL
                      for(i in 1:k)
                      {
                        train <- accepts[ unlist(train_index[i]),]
                        test <- accepts[ unlist(test_index[i]),]#20% test
                        
                        sample_count <-min(summary(as.factor(train$Ac_Y)))
                        print(summary(as.factor(train$Ac_Y)))
                        sample_index <- sample(c(1:length(train$Ac_Y))[which(train$Ac_Y==0)],sample_count)
                        train_DownSample <- data.frame(train[c(sample_index,which(train$Ac_Y==1)),])
                        print(summary(as.factor(train_DownSample$Ac_Y)))
                        print(summary(as.factor(test$Ac_Y)))
                        temp_RF <- randomForest(as.factor(Ac_Y)~.,data = train_DownSample,proximity=TRUE,
                                                na.action =na.omit,mtry=tmp_mtry,ntree=tmp_ntree) #一般对mtry的选择是逐一尝试，直到找到比较理想的值，ntree的选择可通过图形大致判断模型内误差稳定时的值。
                        
                        OOB_err_temp <- c(OOB_err_temp ,mean(temp_RF$err.rate))
                      }
                      
                      OOB_err <- mean(OOB_err_temp)
                      test_mtry <- tmp_mtry
                      test_ntree <-tmp_ntree
                      return(c(OOB_err,test_mtry,test_ntree))
                    }
  stopCluster(cl)

  OOB_err = unlist(result)[seq(1,60,3)]
  test_mtry =  unlist(result)[seq(2,60,3)][which.min(OOB_err)]
  test_ntree = unlist(result)[seq(3,60,3)][which.min(OOB_err)]

  error_matrix <-matrix(0,nrow=20,ncol=8)
  colnames(error_matrix) <- c("MCC","auc","sensitivities","specificities","T0P0","T0P1","T1P0","T1P1")

  for(i in 1:k)
  {
    train <- accepts[ unlist(train_index[i]),]
    test <- accepts[ unlist(test_index[i]),]#20% test
    
    sample_count <-min(summary(as.factor(train$Ac_Y)))
    print(summary(as.factor(train$Ac_Y)))
    sample_index <- sample(c(1:length(train$Ac_Y))[which(train$Ac_Y==0)],sample_count)
    train_DownSample <- data.frame(train[c(sample_index,which(train$Ac_Y==1)),])

    print(summary(as.factor(train_DownSample$Ac_Y)))
    print(summary(as.factor(test$Ac_Y)))

    RF <- randomForest(as.factor(Ac_Y)~.,data = train_DownSample,proximity=TRUE,
                       na.action =na.omit,mtry=test_mtry,ntree=test_ntree)
    save(RF,file = paste0("Result/Model/",cell,"_ActiveP_",Feature_index,"_",i,"RF_fuzz.Rdata"))
    error_matrix[i,1] <- mean(RF$err.rate)
    pred <- predict(RF,test)
    Fred <- t(table(pred,test$Ac_Y))

    modelroc <-roc(test$Ac_Y,as.numeric(pred))
    plot(modelroc,print.auc=TRUE,grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="blue", print.thres=TRUE)

    error_matrix[i,2] <- modelroc$auc
    error_matrix[i,3] <-  modelroc$sensitivities[2]
    error_matrix[i,4] <-  modelroc$specificities[2]
    error_matrix[i,5] <-  Fred[1,1]#Ture0,Pred0  TN
    error_matrix[i,6] <-  Fred[1,2]#Ture0,Pred1  FP
    error_matrix[i,7] <-  Fred[2,1]#Ture1,Pred0  FN
    error_matrix[i,8] <-  Fred[2,2]#Ture1,Pred1  TP
    error_matrix[i,1]<- (Fred[2,2]*Fred[1,1] - Fred[1,2]*Fred[2,1])/(sqrt(Fred[2,2]+Fred[1,2])*sqrt(Fred[2,2]+Fred[2,1])*sqrt(Fred[1,1]+Fred[1,2])*sqrt(Fred[1,1]+Fred[2,1]))
  }
  write.csv(error_matrix,file  = paste0("Result/Summary/",cell,"_Active_P_",Feature_index,"_RF_fuzz.csv"))
}

###########Corss Test RF model 5 Times#############
#'
#' output: test error_matrix saved in "Result/Summary/Coss_"
#'
#' @param model_cell: the cell type name of model
#' @param test_Cell: the cell type name of test
#' @param Feature_index: the feature index of model
#' @param test_Data: the test data about the test cell type
cross_Cell_RF_fuzz <- function(model_cell="GM12878",test_Cell="K562",Feature_index,test_Data)
{
  library("randomForest")
  library(pROC)
  load(test_Data)
  source("Script/make_model_feature_index.R")
  accepts <- Freature_index_make(accepts,Feature_index)

  
  error_matrix <-matrix(0,nrow=20,ncol=4)
  colnames(error_matrix) <- c("","auc","sensitivities","specificities")
  for(i in 1:5)
  {
    load(paste0("Result/Model/",model_cell,"_ActiveP_",Feature_index,"_",i,"RF_fuzz.Rdata"))
    test <- accepts[which(accepts$Ac_Y==0 | accepts$Ac_Y==1),]
    names(error_matrix) <- c("","auc","sensitivities","specificities")

    pred <- predict(RF,test)
    Fred <- t(table(pred,test$Ac_Y))

    modelroc <-roc(test$Ac_Y,as.numeric(pred))
    plot(modelroc,print.auc=TRUE,grid=c(0.1, 0.2),grid.col=c("green", "red"), max.auc.polygon=TRUE,auc.polygon.col="blue", print.thres=TRUE)

    error_matrix[i,2] <- modelroc$auc
    error_matrix[i,3] <-  modelroc$sensitivities[2]
    error_matrix[i,4] <-  modelroc$specificities[2]

  }
  write.csv(error_matrix,file = paste0("Result/Summary/Coss_",model_cell,"_",test_Cell,"_",Feature_index,"_ActiveP_RF_fuzz.csv"))
}





