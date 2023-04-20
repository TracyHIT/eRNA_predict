###################### this function is used to extract the CG% for regions ######################
#' Extract the CG% for regions
#' Output is a Rdata file containing a vector
#'
#' @param regions_list_GRangs a list to get the CG%
#' @param nucfile a file including the GC% of each region
#' @param output_prefix the output file namePrefix containing path
Produce_dataset_pre_nuc <- function(regions_list_GRangs,nucfile="data/MPBS/FANTOM_sorted_delsomeUn_nuc.bed",output_prefix)
{
  print(paste0("Start get the GC_precent in ",length(regions_list_GRangs)," regions!"))
  library("GenomicRanges","ggplot2","eulerr")
  nuc <- read.delim(nucfile, header=FALSE, comment.char="#", stringsAsFactors=FALSE)
  GC_pre_index <- ncol(nuc) - 7
  nuc_gene <- GRanges(seqnames = nuc[,1],ranges = IRanges(start = as.numeric(nuc[,2]),end = as.numeric(nuc[,3])))

  find_match_backgound <- findOverlaps(nuc_gene,regions_list_GRangs)
  GC_precent <- rep(0,length(regions_list_GRangs))
  GC_precent[find_match_backgound@to] <- nuc[find_match_backgound@from,GC_pre_index]
  save(GC_precent,file = paste0(output_prefix,"_nuc_precent.Rdata"))
  print(paste0("End get the GC_precent in ",length(regions_list_GRangs)," regions!"))
}

###################### this function is used to extract the avg_methyl for regions ######################
#' Extract the avg_methyl for regions
#' Output is a Rdata file containing a vector
#'
#' @param regions_list_GRangs a list to get the CG%
#' @param WGBS_file a file including methylation levels on each CpG site
#' @param output_prefix the output file namePrefix containing path
Produce_dataset_avg_region_methyl <- function(regions_list_GRangs,WGBS_file,output_prefix)
{
  print(paste("Start to calculate the avg_methyl for",output_prefix,"."))
  library("GenomicRanges","eulerr")
  library("foreach")
  library("doParallel")
  WGBS <- read.csv(WGBS_file,stringsAsFactors = F)
  WGBS_gene <- GRanges(seqnames = WGBS[,1],ranges = IRanges(start = as.numeric(WGBS[,"position"]),end = as.numeric(WGBS[,"position"])+1))
  WGBS_match <- findOverlaps(regions_list_GRangs,WGBS_gene)

  footprints_index <- tapply(WGBS_match@from, as.factor(WGBS_match@from), function(x){
    return(x[1])
  })
  match_methy<- tapply(WGBS[WGBS_match@to,"methy.read"], as.factor(WGBS_match@from), function(x){
    return(sum(x))
  })
  match_all <- tapply(WGBS[WGBS_match@to,"all.read"], as.factor(WGBS_match@from), function(x){
    return(sum(x))
  })

  methyl <- rep(0,length(regions_list_GRangs))
  methyl[footprints_index] <- match_methy/match_all
  save(methyl,file = paste0(output_prefix,"_avg_methyl.Rdata"))
  print(paste("End to calculate the avg_methyl for",output_prefix,"."))
}

###################### this function is used to extract the L\M\R methylation leves for regions ######################
#' Extract the L\M\R methylation leves for regions
#' Output is a Rdata file containing a matrix ncol=3
#'
#' @param regions_list_GRangs a list to get the CG%
#' @param WGBS_file a file including methylation levels on each CpG site
#' @param output_prefix the output file namePrefix containing path
#' @param regions_list_GRangs_width The region width setted to extract the L\M\R methylation leves.Default is 0，means the ture width of each region in the regions_list_GRangs list.
Produce_dataset_3ad_methyl <- function(regions_list_GRangs,WGBS_file,output_prefix,regions_list_GRangs_width=0)
{
  print(paste("Start to calculate the L_methyl, M_methyl, R_methyl for",output_prefix,"."))
  library("GenomicRanges","ggplot2","eulerr")
  library("foreach")
  library("doParallel")
  WGBS <- read.csv(WGBS_file,stringsAsFactors = F)
  WGBS_gene <- GRanges(seqnames = WGBS[,1],ranges = IRanges(start = as.numeric(WGBS[,"position"]),end = as.numeric(WGBS[,"position"])+1))
  if(regions_list_GRangs_width==0)
  {
    regions_list_GRangs_width <- ceiling(regions_list_GRangs@ranges@width/3)
    R_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = regions_list_GRangs@ranges@start,width =regions_list_GRangs_width ))
    M_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = regions_list_GRangs@ranges@start+regions_list_GRangs_width,width =regions_list_GRangs_width ))
    L_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = regions_list_GRangs@ranges@start+2*regions_list_GRangs_width ,width =regions_list_GRangs_width ))

  }else
  {
    half_width <- ceiling(regions_list_GRangs_width)/2
    midsite <- ceiling(regions_list_GRangs@ranges@width/2)+regions_list_GRangs@ranges@start
    M_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = midsite-half_width,width = regions_list_GRangs_width ))
    L_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = midsite-3*half_width,width = regions_list_GRangs_width ))
    R_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = midsite+half_width,width = regions_list_GRangs_width ))
  }

  WGBS_match_R <- findOverlaps(R_region,WGBS_gene)
  footprints_index <- tapply( WGBS_match_R@from, as.factor( WGBS_match_R@from), function(x){
    return(x[1])
  })
  match_methy<- tapply(WGBS[ WGBS_match_R@to,"methy.read"], as.factor( WGBS_match_R@from), function(x){
    return(sum(x))
  })
  match_all <- tapply(WGBS[ WGBS_match_R@to,"all.read"], as.factor( WGBS_match_R@from), function(x){
    return(sum(x))
  })
  methyl_R <- rep(0,length(regions_list_GRangs))
  methyl_R[footprints_index] <- match_methy/match_all

  WGBS_match_L <- findOverlaps(L_region,WGBS_gene)
  footprints_index <- tapply(WGBS_match_L@from, as.factor( WGBS_match_L@from), function(x){
    return(x[1])
  })
  match_methy<- tapply(WGBS[ WGBS_match_L@to,"methy.read"], as.factor( WGBS_match_L@from), function(x){
    return(sum(x))
  })
  match_all <- tapply(WGBS[WGBS_match_L@to,"all.read"], as.factor(WGBS_match_L@from), function(x){
    return(sum(x))
  })
  methyl_L <- rep(0,length(regions_list_GRangs))
  methyl_L[footprints_index] <- match_methy/match_all

  WGBS_match_M <- findOverlaps(M_region,WGBS_gene)
  footprints_index <- tapply(WGBS_match_M@from, as.factor(WGBS_match_M@from), function(x){
    return(x[1])
  })
  match_methy<- tapply(WGBS[WGBS_match_M@to,"methy.read"], as.factor(WGBS_match_M@from), function(x){
    return(sum(x))
  })
  match_all <- tapply(WGBS[WGBS_match_M@to,"all.read"], as.factor(WGBS_match_M@from), function(x){
    return(sum(x))
  })
  methyl_M <- rep(0,length(regions_list_GRangs))
  methyl_M[footprints_index] <- match_methy/match_all

  region_LMR <- cbind(methyl_L,methyl_M,methyl_R)
  save(region_LMR,file = paste0(output_prefix,"_LMR.Rdata"))
  print(paste("End to calculate the L_methyl, M_methyl, R_methyl for",output_prefix,"."))
}

###################### this function is used to extract the L\M\R histone marker signal score for regions ######################
#' Extract the L\M\R histone marker signals for regions
#' Output is a Rdata file containing a matrix ncol=3
#'
#' @param regions_list_GRangs a list to get the CG%
#' @param bamFile a bam file from ChIP-seq
#' @param output_prefix the output file namePrefix containing path
#' @param regions_list_GRangs_width The region width setted to extract the L\M\R scores.Default is 0，means the ture width of each region in the regions_list_GRangs list.
Produce_dataset_H3K27ac_bam_CPM <-function(regions_list_GRangs,bamFile= "data/Histone_ChIP-seq/ENCFF278LPJ_H3K27ac_HepG2.bam",output_prefix,regions_list_GRangs_width=0)
{
  library(Rsamtools)
  library(rtracklayer)
  library(GenomicAlignments)
  bd <- readGAlignments(bamFile)
  mygr <- as(bd,"GRanges")
  totalReads <- length(mygr)

  if(regions_list_GRangs_width==0)
  {
    regions_list_GRangs_width <- ceiling(regions_list_GRangs@ranges@width/3)
    R_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = regions_list_GRangs@ranges@start,width =regions_list_GRangs_width ))
    M_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = regions_list_GRangs@ranges@start+regions_list_GRangs_width,width =regions_list_GRangs_width ))
    L_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = regions_list_GRangs@ranges@start+2*regions_list_GRangs_width ,width =regions_list_GRangs_width ))

  }else
  {
    half_width <- ceiling(regions_list_GRangs_width)/2
    midsite <- ceiling(regions_list_GRangs@ranges@width/2)+regions_list_GRangs@ranges@start
    M_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = midsite-half_width,width = regions_list_GRangs_width ))
    L_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = midsite-3*half_width,width = regions_list_GRangs_width ))
    R_region <- GRanges(seqnames = regions_list_GRangs@seqnames,ranges = IRanges(start = midsite+half_width,width = regions_list_GRangs_width ))
  }

  bw_match_R <- findOverlaps(R_region,mygr)
  footprints_index <- tapply(bw_match_R@from, as.factor( bw_match_R@from), function(x){
    return(x[1])
  })
  match_score<- tapply(bw_match_R@to, as.factor( bw_match_R@from), function(x){
    return(signif(10^6*length(x)/totalReads,3))
  })
  score_R <- rep(0,length(regions_list_GRangs))
  score_R[footprints_index] <- match_score


  bw_match_L <- findOverlaps(L_region,mygr)
  footprints_index <- tapply( bw_match_L@from, as.factor( bw_match_L@from), function(x){
    return(x[1])
  })
  match_score<- tapply(bw_match_L@to, as.factor(bw_match_L@from), function(x){
    return(signif(10^6*length(x)/totalReads,3))
  })
  score_L <- rep(0,length(regions_list_GRangs))
  score_L[footprints_index] <- match_score



  bw_match_M <- findOverlaps(M_region,mygr)
  footprints_index <- tapply( bw_match_M@from, as.factor( bw_match_M@from), function(x){
    return(x[1])
  })
  match_score<- tapply(bw_match_M@to, as.factor( bw_match_M@from), function(x){
    return(signif(10^6*length(x)/totalReads,3))
  })
  score_M <- rep(0,length(regions_list_GRangs))
  score_M[footprints_index] <- match_score

  score_LMR <- cbind(score_L,score_M,score_R)
  save(score_LMR,file = paste0(output_prefix,"_Score_LMR.Rdata"))
  print(paste("End to calculate the histone LMR for",output_prefix,"."))

}


###################### this function is used to extract the 10 up/downstream geneexpression for regions for regions ######################
#' Extract the 10 up/downstream geneexpression for regions
#' Output is a Rdata file containing a matrix, last colomn is for gene expression
#'
#' @param cell the cell type
#' @param download_gene_file  a .tsv file for gene expression
#' @param gene_Ensemble_annotation a Rdata file for gene Ensemble annotation
get_expression <- function(cell="GM12878",download_gene_file = "data/GM12878_geneExpression_ENCFF912QPL.tsv",
                           gene_Ensemble_annotation="data/gene_Ensemble_annotation.Rdata")
{
  load(gene_Ensemble_annotation)
  GM12878_geneExpression <- read.delim(download_gene_file)
  geneID <- apply(GM12878_geneExpression,1,function(x){strsplit(x[2],split="\\.")[[1]][1]})
  geneExpression <- tapply(GM12878_geneExpression[,"TPM"], as.factor(geneID),FUN = function(x){sum(x)/length(x)})
  geneID <- tapply(geneID,as.factor(geneID),function(x){x[1]})

  cut_gene_num <- 10 #up 10 and down 10
  FANTOM_65407_hg38_delsomeUn <- read.delim("data/MPBS/FANTOM_65407_hg38_delsomeUn.bed", header=FALSE)#65399
  up_down_gene_expression<- apply(FANTOM_65407_hg38_delsomeUn,1,function(x)
  {
    temp <- rep(-1,2*cut_gene_num )
    subgen_anno <- gene_Ensemble_annotation[gene_Ensemble_annotation[,1] == as.character(x[1]),]
    distance <- as.numeric(subgen_anno[,2]) - mean(c(as.numeric(x[2]),as.numeric(x[3])))
    names(distance) <- as.character(subgen_anno[,3])
    distance <- distance[distance < 100000000]
    index <- which.min(abs(distance))
    if(index > cut_gene_num & index < length(distance-cut_gene_num))
    {
      if(distance[index] > 0)
      {
        temp[1:cut_gene_num] <- geneExpression[match(names(distance[(index-cut_gene_num):(index-1)]),geneID)]
        temp[(cut_gene_num+1):(2*cut_gene_num)] <- geneExpression[match(names(distance[(index):(index+cut_gene_num-1)]),geneID)]

      }else
      {
        temp[1:cut_gene_num] <- geneExpression[match(names(distance[(index-cut_gene_num+1):index]),geneID)]
        temp[(cut_gene_num+1):(2*cut_gene_num)] <- geneExpression[match(names(distance[(index+1):(index+cut_gene_num)]),geneID)]
      }
    }
    else
    {print(x)}
    return(temp)
  })
  GM12878_10_geneexpression<- cbind(FANTOM_65407_hg38_delsomeUn,t(up_down_gene_expression))
  save(GM12878_10_geneexpression,file =paste0("data/Feature_X/",cell,"_geneexpression_65407_delsomeUn.Rdata"))
}