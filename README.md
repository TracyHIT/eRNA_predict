# eRNA_predict
# Fitting the data distribution of TPM obtained by CAGE-seq technology
_Script/Fit_noise_From_distribution.R_

# The function for training and testing models
_Script/Ac_enhancer_FANTON_ActiveP_RF_fivefold.R_<br>
_Script/Ac_enhancer_FANTON_ActiveP_XGBoost_fivefold.R_

# Trained models
The performance of RF and XGBoost models is similar. <br> 
Here is the XGBoost model trained on GM12878 data. <br>
Different feature combinations were used. <br>
Please refer to _Script/make_model_feature_index.R_ for the index of feature combinations. <br>
Each feature combination was trained 5 times. Five predictions can be made in parallel, and then integrated through voting. <br>

# data file
FANTOM_Mappd_to_Hg38_genome_loci.Rdata: the gene regions from the FANTOM5, mapped to hg38, 65399 regions. <br>
MPBS/FANTOM_65407_hg38_delsomeUn_nuc.bed: the nucleotide composition for the regions mapped from hg19(65407 regions) to hg38(65399 regions). <br>
gene_Ensemble_annotation.Rdata: gene_Ensemble_annotation,TSS, hg38. <br>
Compare_analysis_data: Store preprocessed data for comparative analysis. <br>
