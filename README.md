# eRNA_predict
# Fitting the data distribution of TPM obtained by CAGE-seq technology
_Script/Fit_noise_From_distribution.R_

# The function for training and testing models
_Script/Ac_enhancer_FANTON_ActiveP_RF_fivefold.R_<br>
_Script/Ac_enhancer_FANTON_ActiveP_XGBoost_fivefold.R_

# Trained models
The performance of RF and XGBoost models is similar. 
Here is the XGBoost model trained on GM12878 data. 
Different feature combinations were used. 
Please refer to _Script/make_model_feature_index.R_ for the index of feature combinations. 
Each feature combination was trained 5 times. Five predictions can be made in parallel, and then integrated through voting.
