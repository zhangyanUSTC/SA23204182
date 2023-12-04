#' @title Binary mixed normal data set
#' @name data_mix2
#' @description Is an example of a data set with a length of 100 and a distribution of 0.6N(0,1)+0.4N(8,4).
#' @examples
#' \dontrun{
#' data(data_mix2)
#' clu<- norm_clus(data_mix2,2)
#' clu$label
#' }
NULL


#' @importFrom Rcpp evalCpp
NULL


#' @title Classify the mixed normal data
#' @description The EM algorithm is used to classify data from binary or ternary positive mixing models.
#' @param data Data set for clustering, length \code{n}
#' @param K EM Indicates the maximum number of iterations of the algorithm
#' @param iters Whether to output EM number of iterations, the default value is T
#' @return A list consisting of the following sections:
#' \describe{
#'   \item{\code{label}}{Data of length \code{n} consisting of 0 or 1 or 2, used to mark the classification category}
#'   \item{\code{param}}{Estimates of the mean and standard deviation of the corresponding two or three normal distributions, as well as estimates of the mixed proportions}
#' }
#' @examples
#' \dontrun{
#' data(data_mix2)
#' clu<- norm_clus(data_mix2,2)
#' sum(clu$label)
#' }
#' @import boot
#' @import MASS
#' @import bootstrap
#' @import DAAG
#' @import coda
#' @import Rcpp
#' @importFrom stats dnorm sd
#' @useDynLib SA23204182
#' @export
norm_clus<-function(data,K=10,iters=T){
  N<-length(data)
  data_bar<-mean(data);data_sd<-sd(data)
  param<-c(0.5,data_bar,data_bar+data_sd,data_sd^2,data_sd^2)
  label<-numeric(N)#label
  for (k in 1:K) {
    for (i in 1:N) {#E
      p1<-param[1]*dnorm(data[i],param[2],sqrt(param[4]))
      p0<-(1-param[1])*dnorm(data[i],param[3],sqrt(param[5]))
      label[i]<-p1/(p1+p0)
    }
    label<-round(label)
    if(k>1 && all(z==label)) {
      if (iters==T) print(c("EM iterations:",k-1))
      break
      }
    z<-label
    param[1]<-sum(label)/N#label=1
    
    param[2]<-crossprod(label,data)/sum(label)
    param[3]<-(sum(data)-crossprod(label,data))/(N-sum(label))
    
    param[4]<-(sum(label*(data-param[2])^2))/sum(label)
    param[5]<-(sum((1-label)*(data-param[3])^2))/(N-sum(label))
    
  }
  z<-round(label)
  return(list(label=z,param=param))
}

#' @title Benchmark R and Rcpp functions.
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of R functions (\code{norm_clus}) and Cpp functions (\code{norm_clus_C}).
#' @examples
#' \dontrun{
#' data(data_mix2)
#' tm1 <- microbenchmark::microbenchmark(
#'   rnR = norm_clus(data_mix2,10,iters=FALSE),
#'   rnC = norm_clusC(data_mix2,10,iters=FALSE))
#' print(summary(tm1)[,c(1,3,5,6)])
#' }
#' @import microbenchmark
NULL

