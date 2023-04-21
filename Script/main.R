source("Script/Make_Traning_Test_dataset.R")
call_generate_X()#To generated the feature data
cell_Make_train_test_data_leftsomeFuzzregionas()#To generated a big matrix data, including all types of features, the training and test data

source("Script/Ac_enhancer_FANTON_ActiveP_RF_fivefold.R")
for(i in c(1:11))#To tain and Test with different feature index within cell type,RF model
{
  for (cell in c("GM12878","K562","HepG2")) {
    data_file <- paste0("Result/Model_ActiveP_data/",cell,"_active_P_data_Filledmiss_fuzz.Rdata")
    print(data_file)
    RF_train_self_tuning_fuzz(data_file,cell,i)
  }
}
for(Feature_index in 1:11)#To test with different feature index corss cell type,RF model
{
  print(paste0("cross ",Feature_index))
  cross_Cell_RF_fuzz("GM12878","K562",Feature_index,paste0("Result/Model_ActiveP_data/K562_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_RF_fuzz("GM12878","HepG2",Feature_index,paste0("Result/Model_ActiveP_data/HepG2_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_RF_fuzz("K562","GM12878",Feature_index,paste0("Result/Model_ActiveP_data/GM12878_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_RF_fuzz("K562","HepG2",Feature_index,paste0("Result/Model_ActiveP_data/HepG2_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_RF_fuzz("HepG2","K562",Feature_index,paste0("Result/Model_ActiveP_data/K562_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_RF_fuzz("HepG2","GM12878",Feature_index,paste0("Result/Model_ActiveP_data/GM12878_active_P_data_Filledmiss_fuzz.Rdata"))
}

source("Script/Ac_enhancer_FANTON_ActiveP_XGBoost_fivefold.R")
for(i in c(1:11))#To tain and Test with different feature index within cell type,XGBoost model
{
  for (cell in c("GM12878","K562","HepG2")) {
    data_file <- paste0("Result/Model_ActiveP_data/",cell,"_active_P_data_Filledmiss_fuzz.Rdata")
    print(data_file)
    xgb_train_self_tuning_fuzz(data_file,cell,i)
  }
}
for(Feature_index in 1:11)#To test with different feature index corss cell type,XGBoost model
{
  print(paste0("cross ",Feature_index))
  cross_Cell_xgb_fuzz("GM12878","K562",Feature_index,paste0("Result/Model_ActiveP_data/K562_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_xgb_fuzz("GM12878","HepG2",Feature_index,paste0("Result/Model_ActiveP_data/HepG2_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_xgb_fuzz("K562","GM12878",Feature_index,paste0("Result/Model_ActiveP_data/GM12878_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_xgb_fuzz("K562","HepG2",Feature_index,paste0("Result/Model_ActiveP_data/HepG2_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_xgb_fuzz("HepG2","K562",Feature_index,paste0("Result/Model_ActiveP_data/K562_active_P_data_Filledmiss_fuzz.Rdata"))
  cross_Cell_xgb_fuzz("HepG2","GM12878",Feature_index,paste0("Result/Model_ActiveP_data/GM12878_active_P_data_Filledmiss_fuzz.Rdata"))
}