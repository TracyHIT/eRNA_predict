###################### this function is used to extract the methyl\histone signal for regions for regions ######################
#' Extract the methyl\histone signal for regions
#' Output is four Rdata file containing avg_methyl,LMR_methyl,LMR_H3K27ac,LMR_H3K9ac
#'
#' @param cell the cell type,here we calculated GM12878,K562,HepG2
#' @param WGBS_file a file including methylation levels on each CpG site
#' @param output_index the output file namePrefix containing path
#' @param regions_file the bedfile saved regions information
#' @param footprint_mpbs_file the bedfile saved motif information
#' @param bamfile_H3K27ac the bamfile saved H3K27ac information
#' @param bamfile_H3K9ac the bamfile saved H3K9ac information
#' @param gene_annotation_file a Rdata file for gene Ensemble annotation
#' @param gene_exression_file  a .tsv file for gene expression
generate_X <- function(cell, WGBS_file,output_index,regions_file="data/MPBS/Random_without_exon_FANTOM_sorted_delsomeUn.bed",
                       bamfile_H3K27ac,bamfile_H3K9ac,regions_list_GRangs_width=0,gene_annotation_file="data/gene_Ensemble_annotation.Rdata",
                       gene_exression_file="data/ENCFF212XWC_K562.tsv")
{
  source('Script/Extracting_input_ features.R')
  library("GenomicRanges")
  regions <- read.delim(regions_file, header=FALSE, stringsAsFactors=FALSE)
  regions_list_GRangs <- GRanges(seqnames = regions[,1],ranges = IRanges(start = as.numeric(regions[,2]),end = as.numeric(regions[,3])))
  rm(regions)
  Produce_dataset_avg_region_methyl(regions_list_GRangs,WGBS_file,output_index)
  Produce_dataset_3ad_methyl(regions_list_GRangs,WGBS_file,output_index,regions_list_GRangs_width)
  Produce_dataset_H3K27ac_bam_CPM(regions_list_GRangs,bamfile_H3K27ac,paste0(output_index,"_H3K27ac"),regions_list_GRangs_width)
  Produce_dataset_H3K27ac_bam_CPM(regions_list_GRangs,bamfile_H3K9ac,paste0(output_index,"_H3K9ac"),regions_list_GRangs_width)
  get_expression(cell,gene_exression_file,gene_annotation_file )
}

call_generate_X <- function()
{
  cell="GM12878"
  WGBS_file <- "data/GM12878.GRCh38.WGBS.csv"
  bamfile_H3K27ac="data/Histone_ChIP-seq/ENCFF948GTC_H3K27ac_GM12878.bam"
  bamfile_H3K9ac="data/Histone_ChIP-seq/ENCFF729HQN_H3K9ac_GM12878.bam"
  output_index <- paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn")
  regions_file <- "data/MPBS/FANTOM_65407_hg38_delsomeUn.bed"
  download_file = "data/GM12878_geneExpression_ENCFF912QPL.tsv"
  annotation="data/gene_Ensemble_annotation.Rdata"
  generate_X(cell, WGBS_file,output_index,regions_file,bamfile_H3K27ac,bamfile_H3K9ac,annotation,download_file)

  cell="HepG2"
  WGBS_file <- "data/HepG2.GRCh38.WGBS.csv"
  bamfile_H3K27ac="data/Histone_ChIP-seq/ENCFF278LPJ_H3K27ac_HepG2.bam"
  bamfile_H3K9ac="data/Histone_ChIP-seq/ENCFF798OBX_H3K9ac_HepG2.bam"
  output_index <- paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn")
  regions_file <- "data/MPBS/FANTOM_65407_hg38_delsomeUn.bed"
  download_file = "data/ENCFF878MHG_HepG2.tsv"
  annotation="data/gene_Ensemble_annotation.Rdata"
  generate_X(cell, WGBS_file,output_index,regions_file,bamfile_H3K27ac,bamfile_H3K9ac,annotation,download_file)

  cell="K562"
  WGBS_file <- "data/K562.GRCh38.WGBS.csv"
  bamfile_H3K27ac = "data/Histone_ChIP-seq/ENCFF600THN_H3K27ac_K562.bam"
  bamfile_H3K9ac ="data/Histone_ChIP-seq/ENCFF763ZGN_H3K9ac_K562.bam"
  output_index <- paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn")
  regions_file <- "data/MPBS/FANTOM_65407_hg38_delsomeUn.bed"
  download_file = "data/ENCFF190NFH_K562.tsv"
  annotation="data/gene_Ensemble_annotation.Rdata"
  generate_X(cell, WGBS_file,output_index,regions_file,bamfile_H3K27ac,bamfile_H3K9ac,annotation,download_file)
}
###############make a big matrix contain all input features AND Y for training and testing for 3 cells
#'
#' make a big matrix contain all input features AND Y for training and testing for 3 cells
#' output is a Rdata file in "Result/Model_ActiveP_data/@cell@_active_P_data_Filledmiss_fuzz.Rdata"
#' @param cell a string for cell name
#' @param test_HMM a bedfile for DnaseI regions
Make_train_test_data_leftsomeFuzzregionas <- function(cell,test_HMM)
{
  source("Script/Fit_noise_From_distribution.R")
  library("GenomicRanges")
  library("ggplot2")
  load(paste0("data/Feature_X/",cell ,"_geneexpression_65407_delsomeUn.Rdata"))
  load(paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn_avg_methyl.Rdata"))
  load(paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn_LMR.Rdata"))
  load(paste0("data/Feature_X/FANTOM_65407_hg38_delsomeUn_nuc_precent.Rdata"))
  load(paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn_Score_LMR.Rdata"))
  GM12878_10_geneexpression_all_feature_data <- cbind(GM12878_10_geneexpression[,-c(1:3)],methyl,region_LMR,GC_precent)
  print(dim(GM12878_10_geneexpression_all_feature_data))

  library(mice)
  imp <-  mice(as.data.frame(lapply(GM12878_10_geneexpression_all_feature_data ,as.numeric)),m=5)
  data1 <- complete(imp,action = 1)
  GM12878_10_geneexpression_all_feature <- cbind(GM12878_10_geneexpression[,c(1:3)],data1)

  ##############make_Y
  GM12878_10_geneexpression_1 <- GM12878_10_geneexpression_all_feature
  Feature_region <- GRanges(seqnames = GM12878_10_geneexpression_1[,1],ranges = IRanges(start = as.numeric(GM12878_10_geneexpression_1[,2]), end = as.numeric(GM12878_10_geneexpression_1[,3])))

  Enhancer <- read.delim(paste0("data/FANTOM/",cell,"_hg38_enhances_phase_1_and_2_expression_tpm_matrix.bed"),header = F)
  tpm <- Enhancer[,c(9,10,11)]
  avg_tpm <- apply(tpm,1, sum)/3
  active_P<- Transfrom_tpm_to_active_p(avg_tpm)
  Enhancer <- Enhancer[which(active_P>0.5),]#####step1

  HMM <- read.delim(test_HMM)
  HMM_region <- GRanges(seqnames = HMM [,1],ranges = IRanges(start = as.numeric(HMM[,2]),end = as.numeric(HMM[,3])))
  caldicate_enhancer_region <- GRanges(seqnames = Enhancer [,1],ranges = IRanges(start = as.numeric(Enhancer[,2]),end = as.numeric(Enhancer[,3])))
  HMM_overlap <- findOverlaps(HMM_region,caldicate_enhancer_region)
  Ac_Enhancer_region <- caldicate_enhancer_region[unique(HMM_overlap@to)]

  index <- findOverlaps(Ac_Enhancer_region,Feature_region)#####step2
  Ac_Y <- rep(0,length(Feature_region))
  Ac_Y[index@to] <- 1

  index_0 <- findOverlaps(caldicate_enhancer_region,Feature_region)#####step2,2,3,4
  index_2 <- findOverlaps(HMM_region,Feature_region)
  Ac_Y[which(Ac_Y == 0)[which(Ac_Y == 0) %in% index_2@to]] <- -1
  Ac_Y[which(Ac_Y == 0)[which(Ac_Y == 0) %in% index_0@to]] <- -1


  load(paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn_H3K27ac_bw_Score_LMR.Rdata"))
  H3K27ac_score_LMR_bw<- score_LMR
  load(paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn_H3K27ac_Score_LMR.Rdata"))
  H3K27ac_score_LMR<- score_LMR
  load(paste0("data/Feature_X/",cell,"_FANTOM_65407_hg38_delsomeUn_H3K9ac_Score_LMR.Rdata"))
  H3K9ac_score_LMR<- score_LMR
  #GM12878_10_geneexpression[,-c(1:3)],methyl,region_LMR,GC_precent
  accepts <- data.frame(cbind(GM12878_10_geneexpression_all_feature[,-c(1:3)],H3K27ac_score_LMR_bw,H3K27ac_score_LMR,H3K9ac_score_LMR,Ac_Y))

  save(accepts,file = paste0("Result/Model_ActiveP_data/",cell,"_active_P_data_Filledmiss_fuzz.Rdata"))
}



cell_Make_train_test_data_leftsomeFuzzregionas<-function()
{
  Make_train_test_data_leftsomeFuzzregionas ("GM12878","data/DNaseIAndHmm/wgEncodeAwgDnaseUwdukeGm12878.bed")
  Make_train_test_data_leftsomeFuzzregionas ("K562","data/DNaseIAndHmm/wgEncodeAwgDnaseUwdukeK562.bed")
  Make_train_test_data_leftsomeFuzzregionas ("HepG2","data/DNaseIAndHmm/wgEncodeAwgDnaseUwdukeHepg2.bed")
}


