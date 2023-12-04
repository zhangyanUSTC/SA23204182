#include <Rcpp.h>
using namespace Rcpp;

//' @title Classify the mixed normal data using Rcpp
//' @description The EM algorithm is used to classify data from binary or ternary positive mixing models Rcpp
//' @param data Data set for clustering, length \code{n}
//' @param K EM Indicates the maximum number of iterations of the algorithm
//' @param iters Whether to output EM number of iterations, the default value is T
//' @return A list consisting of the following sections:
//' \describe{
//'   \item{\code{label}}{Data of length \code{n} consisting of 0 or 1 or 2, used to mark the classification category}
//'   \item{\code{param}}{Estimates of the mean and standard deviation of the corresponding two or three normal distributions, as well as estimates of the mixed proportions}
//' }
//' @examples
//' \dontrun{
//' data(data_mix2)
//' norm_clusC(data_mix2,10)
//' }
//' 
//' @export
// [[Rcpp::export]]
List norm_clusC(NumericVector data,int K,bool iters=true) {
  int N=data.size();
  double data_mean=mean(data);
  double data_sd=sd(data);
  NumericVector label(N);
  NumericVector param(5);
  param[0] = 0.5;
  param[1] = data_mean;
  param[2] = data_mean + data_sd;
  param[3] = pow(data_sd, 2);
  param[4] = pow(data_sd, 2);
  NumericVector z(N);
  for (int iter = 0; iter < K; iter++) {
    double p1, p0;
    for (int i = 0; i < N; i++) {
      p1 = param[0]*exp(-0.5 * pow((data[i] - param[1]), 2)/param[3])/ sqrt(param[3]);
      p0 = (1-param[0])*exp(-0.5 * pow((data[i] - param[2]), 2)/ param[4]) / sqrt(param[4]);
      label[i]=p1/(p1+p0);
    }
    for (int i = 0; i < N; i++){
      label[i]=round(label[i]);
    }
    if (iter >0 && is_true(all(z==label))){
      if (iters==true ) 
        Rcout << "EM iterations:" << iter << std::endl; 
      break;
    }
    for (int i = 0; i < N; i++){
      z[i]=label[i];
    }
    param[0]=mean(label);
    double s1=0,s2=0,s3=0,s4=0;
    
    for (int i = 0; i < N; i++) {
      s1 += label[i] * data[i]/N/param[0];
      s2 += (1-label[i]) * data[i]/N/(1-param[0]);
      s3 += label[i] * pow(data[i],2);
      s4 += (1-label[i]) * pow(data[i],2);
    }
    param[1] = s1;
    param[2] = s2;
    param[3] = s3/N/(param[0])- pow(s1,2);
    param[4] = s4/N/(1-param[0])- pow(s2,2);
  }
  List result = List::create(
    Named("label") = label,
    Named("param") = param);
  return result;
}



