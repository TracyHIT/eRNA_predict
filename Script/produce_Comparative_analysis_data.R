#For Enhancer-LSTMAtt, get the sequences based on the test data, length=200bp, center
#first get all the seq
Enhancer_LSTMAtt_input_65399_region_seq <- function(Mappd_to_Hg38_genome_loci_file="data/FANTOM_Mappd_to_Hg38_genome_loci.Rdata",dir = "/home/zqluoximei/ref_hg38/hg38/",flank_size = 100)
{
  library("Biostrings")
  library("rlist")
  library("foreach")
  library("doParallel")
  load(Mappd_to_Hg38_genome_loci_file)
  center_chr_loci <- cbind(Mappd_to_Hg38_genome_loci[,1],
                          ceiling((Mappd_to_Hg38_genome_loci[,2]+Mappd_to_Hg38_genome_loci[,3])/2))
  chr <- ""
  s <- NULL
  Seq <- rep("", nrow(center_chr_loci))
  cl <- makeCluster(8)
  registerDoParallel(cl)
  DNAtens <- foreach(j = c(1:nrow(center_chr_loci))) %dopar% {
    library("Biostrings")
    if (center_chr_loci[j, 1] != chr) {
      chr <- as.character(center_chr_loci[j, 1])
      s <- readDNAStringSet(paste0(dir, chr, ".fa"))
    }
    Seq[j] <- toString(subseq(s, start = as.numeric(as.character(center_chr_loci[j, 2])) - flank_size, end = as.numeric(as.character(center_chr_loci[j, 2])) + flank_size-1))
    return(Seq[j])
  }
  stopCluster(cl)
  Seq <-as.character(unlist(DNAtens))
  all_Seq <- cbind(Mappd_to_Hg38_genome_loci,Seq)
  print(length(Seq))
  colnames(all_Seq) <- c("Chr","start","end","Seq")
  save(all_Seq,file = "data/Compare_analysis_data/Enhancer_LSTMAtt_200_all_seq.Rdata")
}
#Enhancer_LSTMAtt_input_65399_region_seq()

#for each cell,we select 20% regions to make Test data, test 5 times to get the averge ACC Sp ... ...
Enhancer_LSTMAtt_cell_test <- function(cell="GM12878",all_Seq="data/Compare_analysis_data/Enhancer_LSTMAtt_200_all_seq.Rdata")
{
  load(all_Seq)
  load(paste0("Result/Model_ActiveP_data/",cell,"_active_P_data_Filledmiss_fuzz.Rdata"))
  k <- 5
  len <- nrow(all_Seq)
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
  test_index =  t_mysplit$two

  for(i in 1:k)
  {
    temp_index <- unlist(test_index[i])
    test_Y <- accepts$Ac_Y[ temp_index]
    test_X_seq <- all_Seq[temp_index,"Seq"]
    test_X_name <- paste0(">",all_Seq[temp_index,1],"_",all_Seq[temp_index,2],"_",all_Seq[temp_index,3])
    Postive_seq <- as.character(test_X_seq[which(test_Y==1)])
    Postive_name <- test_X_name[which(test_Y==1)]
    Postive <- rep("",2*length(Postive_seq))
    Postive[seq(1,2*length(Postive_seq),2)] <- Postive_name
    Postive[seq(2,2*length(Postive_seq),2)] <- Postive_seq
    write.table(Postive,file = paste0("data/Compare_analysis_data/Enhancer_LSTMAtt_",cell,"_Postive_",i,".txt"),row.names = F,col.names = F,quote = F)

    Negtive_seq <- as.character(test_X_seq[which(test_Y==0)])
    Negtive_name <- test_X_name[which(test_Y==0)]
    Negtive <- rep("",2*length(Negtive_seq))
    Negtive[seq(1,2*length(Negtive_seq),2)] <- Negtive_name
    Negtive[seq(2,2*length(Negtive_seq),2)] <- Negtive_seq
    write.table(Negtive,file = paste0("data/Compare_analysis_data/Enhancer_LSTMAtt_",cell,"_Negtive_",i,".txt"),row.names = F,col.names = F,quote = F)
  }

}
Enhancer_LSTMAtt_cell_test("GM12878")
Enhancer_LSTMAtt_cell_test("K562")
Enhancer_LSTMAtt_cell_test("HepG2")


