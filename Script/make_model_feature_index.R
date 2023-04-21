######select the feature based on the index from the Big matrix including all types of features
#' This is an embedded function that supports model training an
#' output is the feature matrix for training and test based on the feature index
#'
#' @param accepts the big matrix including all types of features
#' @param Feature_index the index of the feature
Freature_index_make <- function(accepts,Feature_index)
{
  max_gene <- apply(accepts[,c(1:20)],1, max)
  max_gene[max_gene==-1]<-0
  max_gene[max_gene>100]<-100
  accepts <- cbind(max_gene,accepts)
  #max_gene,1:20,avg_methyl methyLMR*3,GC_precent,bigwig*3,H3K27ac*3,H3K9ac*3,Ac_Y
  if(Feature_index == 1)#WGBS,GC_precent
  {
    accepts <- accepts[,c(22:26,36)]
  }else
  {
    if(Feature_index == 2)#max_gene,1:20,avg_methyl methyLMR*3,GC_precent
    {
      accepts <- accepts[,c(1:26,36)]
    }else{
      if(Feature_index == 3)#GC_precent,H3K27ac
      {
        accepts <- accepts[,c(26,30:32,36)]
      }else{if(Feature_index == 4)#GC_precent,H3K9ac
        {
          accepts <- accepts[,c(26,33:35,36)]
        }else{if(Feature_index == 5)#GC_precent,H3K27ac,H3K9ac
          {
            accepts <- accepts[,c(26,30:35,36)]
          }else{if(Feature_index == 6)#WGBS,GC_precent,H3K27ac
            {
              accepts <- accepts[,c(22:26,30:32,36)]
            }else {if(Feature_index == 7)#WGBS,GC_precent,H3K9ac
              {
                accepts <- accepts[,c(22:26,33:35,36)]
              }
              else{if(Feature_index == 8)#WGBS,GC_precent,H3K27a,cH3K9ac
                {
                  accepts <- accepts[,c(22:26,30:35,36)]
                }else{if(Feature_index == 9)#gene,WGBS,GC_precent,H3K27ac
                  {
                    accepts <- accepts[,c(1:26,30:32,36)]
                  }else{if(Feature_index == 10)#WGBS,gene,GC_precent,H3K9ac
                    {
                      accepts <- accepts[,c(1:26,33:35,36)]
                    }else{if(Feature_index == 11)#WGBS,gene,GC_precent,H3K9ac,H3K27ac
                      {
                        accepts <- accepts[,c(1:26,30:35,36)]
                      }
                    }
                    
                    
                  }
                }
                
              }
            }
          }
        }
        
      }
    }
  }
  return(accepts)
}