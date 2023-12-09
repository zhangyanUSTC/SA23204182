## ----eval=TRUE----------------------------------------------------------------
library(SA23204182)
data(data_mix2)
clu<-norm_clus(data=data_mix2 ,K=10)

## -----------------------------------------------------------------------------
str(clu)
clu$label

## ----eval=FALSE---------------------------------------------------------------
#  norm_clus<-function(data,K =10){
#    N<-length(data)
#    data_bar<-mean(data);data_sd<-sd(data)
#    param<-c(0.5,data_bar,data_bar+data_sd^2/2,data_sd^2,data_sd^2)
#    label<-numeric(N)#label
#    for (k in 1:K) {
#      for (i in 1:N) {#E
#        p1<-param[1]*dnorm(data[i],param[2],sqrt(param[4]))
#        p0<-(1-param[1])*dnorm(data[i],param[3],sqrt(param[5]))
#        label[i]<-p1/(p1+p0)
#      }
#      label<-round(label)
#      if(k>1 && all(z==label)) {print(c("EM iterations:" ));break}
#      z<-label
#      param[1]<-sum(label)/N#label=1
#  
#      param[2]<-crossprod(label,data)/sum(label)
#      param[3]<-(sum(data)-crossprod(label,data))/(N-sum(label))
#  
#      param[4]<-(sum(label*(data-param[2])^2))/sum(label)
#      param[5]<-(sum((1-label)*(data-param[3])^2))/(N-sum(label))
#  
#    }
#    z<-round(label)
#    return(list(label=z,param=param))
#  }

## ----eval=TRUE----------------------------------------------------------------
# library(SA23204182)
data(data_mix2)
cluC<-norm_clusC(data=data_mix2 ,K=10)

## -----------------------------------------------------------------------------
str(cluC)
cluC$label

## ----eval=FALSE---------------------------------------------------------------
#  List norm_clusC(NumericVector data,int K) {
#    int N=data.size();
#    double data_mean=mean(data);
#    double data_sd=sd(data);
#    NumericVector label(N);
#    NumericVector param(5);
#    param[0] = 0.5;
#    param[1] = data_mean;
#    param[2] = data_mean + data_sd;
#    param[3] = pow(data_sd, 2);
#    param[4] = pow(data_sd, 2);
#    NumericVector z(N);
#    for (int iter = 0; iter < K; iter++) {
#      double p1, p0;
#      for (int i = 0; i < N; i++) {
#        p1 = param[0]*exp(-0.5 * pow((data[i] - param[1]), 2)/param[3])/ sqrt(param[3]);
#        p0 = (1-param[0])*exp(-0.5 * pow((data[i] - param[2]), 2)/ param[4]) / sqrt(param[4]);
#        label[i]=p1/(p1+p0);
#      }
#      for (int i = 0; i < N; i++){
#        label[i]=round(label[i]);
#      }
#      if (iter >0 && is_true(all(z==label))){
#        Rcout << "EM iterations:" << iter << std::endl;
#        break;
#      }
#      for (int i = 0; i < N; i++){
#        z[i]=label[i];
#      }
#      param[0]=mean(label);
#      double s1=0,s2=0,s3=0,s4=0;
#  
#      for (int i = 0; i < N; i++) {
#        s1 += label[i] * data[i]/N/param[0];
#        s2 += (1-label[i]) * data[i]/N/(1-param[0]);
#        s3 += label[i] * pow(data[i],2);
#        s4 += (1-label[i]) * pow(data[i],2);
#      }
#      param[1] = s1;
#      param[2] = s2;
#      param[3] = s3/N/(param[0])- pow(s1,2);
#      param[4] = s4/N/(1-param[0])- pow(s2,2);
#    }
#    List result = List::create(
#      Named("label") = label,
#      Named("param") = param);
#    return result;
#  }

## ----eval=T-------------------------------------------------------------------
# library(SA23204182)
data(data_mix2)
tm1 <- microbenchmark::microbenchmark(
rnR = norm_clus(data_mix2,10,iters=FALSE),
rnC = norm_clusC(data_mix2,10,iters=FALSE))
print(summary(tm1)[,c(1,3,5,6)])

## ----eval=FALSE---------------------------------------------------------------
#  n <- 100
#  adjacency_Matrix <- matrix(rbinom(n^2, 1, 0.1), n, n)

## -----------------------------------------------------------------------------
data(adjacency_matrix)

## -----------------------------------------------------------------------------
numClusters <- 2
result <- spectralClustering(adjacencyMatrix, numClusters)
cat("Community Assignments:", result, "\n")

## -----------------------------------------------------------------------------
library(igraph)
g <- make_graph('Zachary')
plot(g, layout = layout_with_kk)

## -----------------------------------------------------------------------------
adjacencyMatrixg <- get.adjacency(g)
numClusters <- 2
result <- spectralClustering(adjacencyMatrixg, numClusters)
cat("Community Assignments:", result, "\n")

## ----eval=F-------------------------------------------------------------------
#  spectralClustering <- function(adjacencyMatrix, numClusters) {
#    graph <- graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected", weighted = NULL)
#  
#    # Get the Laplacian matrix
#    laplacianMatrix <- graph.laplacian(graph, normalized = TRUE)
#  
#    # Compute the eigenvalues and eigenvectors of the Laplacian matrix
#    laplacianEigen <- eigen(laplacianMatrix)
#  
#    # Select the first numClusters eigenvectors
#    selectedEigenVectors <- laplacianEigen$vectors[, 1:numClusters]
#  
#    # Perform spectral clustering
#    spectralClusters <- kmeans(selectedEigenVectors, centers = numClusters)$cluster
#  
#    return(spectralClusters)
#  }

