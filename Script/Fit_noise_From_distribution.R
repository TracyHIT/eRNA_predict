###########Transfrom_tpm_to_active_p
#'
#' Noise obeys Poisson distribution， true eRNA signal debinormal Distribution
#' output is a vector of eRNA probability. The length of the vector is the same as the length of the input vector
#'
#' @param avg_tpm average tpm of three replicates
#' @param figure_index the prefix of figure,default is "",mains do not save figure
Transfrom_tpm_to_active_p <- function(avg_tpm,figure_index="")
{
  active_p <- rep(0,length(avg_tpm))
  active_p[avg_tpm<0.001]=0
  num <- ceiling(10*(logb(avg_tpm[avg_tpm>0.001])+4))#########we have to get the avg_tpm to this

  hist(num,freq = F)

  library(MASS)
  l1 = function(para)
  {
    f1 = dnbinom(num,para[1],1/(1+para[2]))
    f2 = dpois(num,para[3],log = FALSE)
    f= para[4]*f1+(1-para[4])*f2
    l1 = sum(log(f))
    return(-l1)
  }
  geyser.est=nlminb(c(4,6,4,0.4),l1,lower=c(0.0001,-Inf,0.0001,0.3),upper=c(200,Inf,Inf,0.9))
  p <- geyser.est$par
  p
  x <- seq(0,80)
  f = p[4]*dnbinom(x,p[1],1/(1+p[2]))+(1-p[4])*dpois(x,p[3],log = FALSE)
  ############plot

  if(figure_index=="")
  {}else
  {
    pdf(file = paste0(figure_index,"_distribution_active_p.pdf"))
    hist(num,freq = F)
    lines(x,f)
    dev.off()}

  num <- ceiling(10*(logb(avg_tpm[avg_tpm>0.001])+4))
  PP<- dnbinom(num,p[1],1/(1+p[2]))
  PN <- dpois(num,p[3],log = FALSE)
  active_p  <- PP/(PP+PN)
  active_p[which(is.na(active_p))] <- 0
  return(active_p)
}