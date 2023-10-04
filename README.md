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

# data file
_FANTOM_Mappd_to_Hg38_genome_loci.Rdata: the gene regions from the FANTOM5, mapped to hg38, 65399 regions.
_MPBS/FANTOM_65407_hg38_delsomeUn_nuc.bed: the nucleotide composition for the regions mapped from hg19(65407 regions) to hg38(65399 regions).
_gene_Ensemble_annotation.Rdata: gene_Ensemble_annotation,TSS, hg38.
_Compare_analysis_data: Store preprocessed data for comparative analysis.
