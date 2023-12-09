## -----------------------------------------------------------------------------
set.seed(111)
indept.Metropolis <- function(sigma_g, x0, N){#链生成
  #sigma_g为提议方差 x0为初始值 N链长 
    x = numeric(N)
    x[1] = x0
    u = runif(N)
    k = 0
    for(i in 2:N){
      y = rnorm(1, x[i-1], sigma_g)
      porb = (1+x[i-1]^2) / (1+y^2)
      if(u[i] <= porb)
        x[i] = y  
      else{
        x[i] = x[i-1]
        k = k + 1
      }
    }
    return(list(x=x, k=k))
  }
chain_len = 5000
xt = numeric(chain_len)
sigma_g = c(0.05, 0.5, 2.5, 16)
x0 <- 3   # initial point
  indept1 <- indept.Metropolis(sigma_g[1], x0, chain_len)
  indept2 <- indept.Metropolis(sigma_g[2], x0, chain_len)
  indept3 <- indept.Metropolis(sigma_g[3], x0, chain_len)
  indept4 <- indept.Metropolis(sigma_g[4], x0, chain_len)
print(c(indept1$k, indept2$k, indept3$k, indept4$k)/chain_len)#四个链的拒绝率
  refline <- qcauchy(c(.025, .975))
  indept <- cbind(indept1$x, indept2$x, indept3$x, indept4$x)
  for (j in 1:4) {
      plot(indept[,j], type="l",
           xlab=bquote(sigma == .(round(sigma_g[j],3))),
           ylab="X", ylim=range(indept[,j]))
      abline(h=refline)
  }

## -----------------------------------------------------------------------------
ffpi<-function(theta){
  1/(1+theta^2)*exp(-(theta-6)^2/2)
}
piii<-function(theta){
  z <- integrate(ffpi,-Inf,Inf)$value
  y <- ffpi(theta)/z
  y
}#后验分布
curve(piii,-5,10)
#单峰

## -----------------------------------------------------------------------------
quantt1 <- function(x){
   integrate(piii,-Inf,x)$value-0.025
}
quantt2 <- function(x){
   integrate(piii,-Inf,x)$value-0.975
}
x025<-uniroot(quantt1,c(-5,10))$root
x975<-uniroot(quantt2,c(-5,10))$root
qq1<-piii(x025)
qq2<-piii(x975)
x025;x975

## -----------------------------------------------------------------------------
repeat{
kk<-(qq1+qq2)/2
ppkk <- function(theta,kk){
  y <- piii(theta)-kk
  y
}
cc1<- uniroot(ppkk,c(0,6),kk=kk)$root
cc2<- uniroot(ppkk,c(6,10),kk=kk)$root
ppp<- integrate(piii,cc1,cc2)$value
   if (abs(ppp-0.95)<=0.0001) break
   if (p-0.95>0)  qq2<-kk
   else  qq1<-kk
   } 
ppp  #HPD可信区间的可信度

cc1;cc2 

## -----------------------------------------------------------------------------
my.sample1 =function(x,pro,size){#总体 概率 样本容量
  cdfpro<- cumsum(pro); m <- size; U = runif(m)#根据pro求分布函数Fx,并生成U
  r <- x[findInterval(U,cdfpro)+1]#逆变换
  ct <- as.vector(table(r))
  print(c("各元素被抽中频率：",ct/m))
  return(r)
}

x <- c("a","b","c","d"); pro <- c(.1,.1,.3,.5)#抽样总体与概率
my.sample1(x,pro,100)

## -----------------------------------------------------------------------------
my.sample2 = function(x,size,Lower_bound){#条件抽样 抽高于Lower_bound的元素
  xr<-0
  xx<-c()#x中符合条件的角标
  for (i in 1:length(x)){
    if (x[i]>Lower_bound){
      xr=xr+1
      xx<-c(xx,i)
    }
  }
  if (length(xx)==0)
    stop("没有符合条件的元素！")
  U = runif(size)#抽样数目
  cdfp=1:length(xx)/length(xx)
  return(x[xx[findInterval(U,cdfp)+1]])
}

x <- c(1,3,5,7,9,8)#抽样总体
my.sample2(x,100,5)

## -----------------------------------------------------------------------------
fni=function(x){
  if (x<=0.5)
    return(log(2*x))
  if (x>0.5)
    return(-log(2*(1-x)))
}#逆函数
x<-c()
U = runif(1000)
for (i in U) {
  x<-c(x,fni(i))
}#代入逆函数
minx =floor(min(x))
maxx = floor(max(x))+1
hist(x,probability = TRUE,ylim = c(0,0.6),xlim = c(minx,maxx),
     breaks = seq(minx,maxx,by = 0.2),col = "skyblue")
y <- seq(-10,10,0.1)
lines(y,0.5*exp(-abs(y)),lwd = 2) #绘制直方图与密度函数图像进行对比

## -----------------------------------------------------------------------------
n <- 1e3;j<-k<-0;y <- numeric(n)#k样本数，y样本
while (k < n) {
  u <- runif(1)
  j <- j + 1#试验次数
  x <- runif(1) #产生新候选样本x
  if (12*x*x*(1-x)>2*u) {
    #接受x
    k <- k + 1
    y[k] <- x
  }
}
hist(y,probability = TRUE,ylim = c(0,3),xlim = c(0,1),breaks = seq(0,1,by = 0.04),col = "skyblue")
z <- seq(0,1,0.01)
lines(z,12*z*z*(1-z),lwd = 2,col="red") #绘制直方图与pdf曲线
print(c("接受率为：",k/j))

## -----------------------------------------------------------------------------
n <- 1e4;j<-k<-0;x <- numeric(n)#k样本数，x样本
while (k < n) {
  u<- runif(3,min=-1,max=1)
  j <- j + 1#试验次数
  if (abs(u[3])>abs(u[2]) & abs(u[3])>abs(u[1])) {
    #接受u2
    k <- k + 1
    x[k] <- u[2]
  }
  else {
    k <- k + 1
    x[k] <- u[3]
  }
}
density_e <- density(x)
hist(x,probability = TRUE,ylim = c(0,1),xlim = c(-1,1),breaks = seq(-1,1,by = 0.04),col = "skyblue")
y <- seq(-1,1,0.01)
lines(y,0.75*(1-y*y),lwd = 2) #绘制直方图与pdf曲线
lines(density_e, col = "red", lwd = 2)
legend("topright", legend = "Estimate", col = "red", lwd = 2)
legend("topleft", legend = "pdf", col = "black", lwd = 2)

## -----------------------------------------------------------------------------
set.seed(1234)
l <- c(0.5,0.8,1)
d <- 1
m <- 1e6
k<-100
X <- array(runif(m*k,0,d/2),dim=c(m,k))#随机数序号对应行、实验序号对应不同列
Y <- array(runif(m*k,0,pi/2),dim=c(m,k))
some<-outer(l/2,sin(Y))>outer(c(1,1,1),X)#张量积，在前面补充参数维数
#apply(some,c(1,3),mean)#1不同参数.2单次实验随机数。3.重复实验次数
pihat<-2*outer(l,array(1,dim = k))/apply(some,c(1,3),mean)#求估计
varpi<-apply(pihat,1,var)
varpi

## -----------------------------------------------------------------------------
varth<-(-1)*(pi)^2*(2*l-pi)/2/m/l
varth

## -----------------------------------------------------------------------------
abs(varth-varpi)/varth

## -----------------------------------------------------------------------------
set.seed(1234)
m <- 1000
k=100
x <- array(runif(m*k),dim=c(m,k))
theta_anti <- apply((exp(x)+exp(1-x)) / 2,2,mean)#theta_anti <- apply((exp(x)+exp(1-x))[1:m/2,1:k]/ 2,2,mean)
theta_sim <- apply(exp(x),2,mean)
(var(theta_sim)-var(theta_anti))/var(theta_sim)



## -----------------------------------------------------------------------------
#图像比对
x <- seq(1, 5, by = 0.1)
y1 <- x^2/sqrt(2*pi)*exp(-x^2/2)
y2 <- 2*exp(-(x-1)^2/2)/sqrt(2*pi)
y3 <- 1/x^2
  
plot(x, y1, type = "l", col = "blue", xlab = "x", ylab = "y",ylim=c(0,1))
lines(x, y2, col = "red")
lines(x, y3, col = "green")
legend("topright", legend = c("g(x)", "f1(x)","f2(x)"), col = c("blue", "red","green"), lty = 1)

## -----------------------------------------------------------------------------
f1<-function(x){
  return(x^4/2/sqrt(2*pi)*exp(-x^2/2-x+0.5))}
result1 <- integrate(f1, lower = 1, upper = Inf)
result1#f1对应的积分

f2<-function(x){
  return(x^6/(2*pi)*exp(-x^2))}
result2 <- integrate(f2, lower = 1, upper = Inf)
result2#f2对应的积分

## -----------------------------------------------------------------------------
set.seed(123)
x<-rnorm(100)
y<-abs(x-1)+1
gf<-function(x){#函数g/f
  return(x^2/sqrt(2*pi)*exp(-x^2/2)/(2*exp(-(x-1)^2/2)/sqrt(2*pi)))}
theta.hat<-mean(gf(y))

g<-function(x){#函数g(x)
  return( x^2/sqrt(2*pi)*exp(-x^2/2))
}
theta.real<-integrate(g, lower = 1, upper = Inf)
theta.hat;theta.real

## -----------------------------------------------------------------------------
set.seed(1)
M <- 10000
g <- function(x) {#函数g(x)
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}
u <- c(runif(M/5,0,0.2),runif(M/5,0.2,0.4),runif(M/5,0.4,0.6),runif(M/5,0.6,0.8),runif(M/5,0.8,1.0))#分层抽样

x<--log(1-u*(1-exp(-1)))
fg <- g(x)/exp(-x)*(1-exp(-1))

theta.hat <- mean(fg)#重要函数估计
fgr<-apply(array(fg,dim=c(M/5,5)),1,sum)/5#分层组合
se <- sd(fgr)
theta.hat;se


## -----------------------------------------------------------------------------
set.seed(123)
m=10000
n <- 20
alpha <- .05
torf<-numeric(m)
for (i in 1:m) {
  x <-rchisq(n, df=2)#rnorm(n)+2#
  torf[i] <-abs(2-mean(x)) <qt(1-alpha/2,df=n-1)*sqrt(var(x))/sqrt(n)#是否落入置信区间
}
sum(torf)/m

## -----------------------------------------------------------------------------
n<-20;m<-10000;p1<-p2<-p3<-numeric(m)
for (i in 1:m) {#第一种卡方分布
  x<-rchisq(n,df=1)
  pz<-t.test(x,mu=1)
  p1[i]<-pz$p.value
}
mean(p1<0.05)#第一类错误率

for (i in 1:m) {#第二种均匀分布
  x<-2*runif(n)
  pz<-t.test(x,mu=1)
  p2[i]<-pz$p.value
}
mean(p2<0.05)#第一类错误率

for (i in 1:m) {#第三种指数分布
  x<-rexp(n)
  pz<-t.test(x,mu=1)
  p3[i]<-pz$p.value
}
mean(p3<0.05)#第一类错误率


## -----------------------------------------------------------------------------
n<-100;m<-10000;p1<-p2<-p3<-numeric(m)
for (i in 1:m) {#第一种卡方分布
  x<-rchisq(n,df=1)
  pz<-t.test(x,mu=1)
  p1[i]<-pz$p.value
}
mean(p1<0.05)#第一类错误率

for (i in 1:m) {#第二种均匀分布
  x<-2*runif(n)
  pz<-t.test(x,mu=1)
  p2[i]<-pz$p.value
}
mean(p2<0.05)#第一类错误率

for (i in 1:m) {#第三种指数分布
  x<-rexp(n)
  pz<-t.test(x,mu=1)
  p3[i]<-pz$p.value
}
mean(p3<0.05)#第一类错误率


## -----------------------------------------------------------------------------
#setting
set.seed(1)
m<-1000;M<-1000
FWER1<-numeric(M);FDR1<-numeric(M);TPR1<-numeric(M)
FWER2<-numeric(M);FDR2<-numeric(M);TPR2<-numeric(M)
#compute
p<-numeric(m)#记录p值
for (i in 1:M) {
  p[1:950]<-runif(m*.95);p[951:m]<-rbeta(m*.05,0.1,1)
  p.adj1 = p.adjust(p,method='bonferroni')#p*N
  p.adj2 = p.adjust(p,method='fdr')#p*N/max(k)(就是有相同值时取大秩)
  #计算各值
  FWER1[i]<-sum(p.adj1[1:m*.95]<0.1)>0
  FWER2[i]<-sum(p.adj2[1:m*.95]<0.1)>0
  FDR1[i]<-sum(p.adj1[1:m*.95]<0.1)/sum(p.adj1<0.1)
  FDR2[i]<-sum(p.adj2[1:m*.95]<0.1)/sum(p.adj2<0.1)
  TPR1[i]<-sum(p.adj1[951:m]<0.1)/50
  TPR2[i]<-sum(p.adj2[951:m]<0.1)/50
}

## ----echo=FALSE---------------------------------------------------------------
A <- matrix(round(c(sum(FWER1)/M,sum(FWER2)/M,mean(FDR1),mean(FDR2),mean(TPR1),mean(TPR2)),3),2)
colnames(A) <- c('FWER','FDR','TPR')
rownames(A) <- c('Bonferroni矫正','B-H矫正')
knitr::kable(A, format = "markdown")

## ----eval=FALSE---------------------------------------------------------------
#  #setting
#  set.seed(1)
#  lammda<-2
#  n<-c(5,10,20);B<-1e3;m<-1e3
#  bias<-matrix(0, nrow =m, ncol = 3)
#  se<-matrix(0, nrow =m, ncol = 3)
#  
#  #n=5时
#  for (k in 1:m) {
#    x <- rexp(n[1],lammda)
#    lamstar<-numeric(B)
#    for(b in 1:B){
#      xstar <- sample(x,replace=TRUE)
#      lamstar[b]<-1/mean(xstar)
#    }
#    bias[k,1]<-mean(lamstar)-1/mean(x)
#    se[k,1]<-sd(lamstar)
#  }
#  round(c(bias1=mean(bias[,1]),se1=mean(se[,1]),bias1.theo=lammda/(n[1]-1),se1.theo=lammda*n[1]/(n[1]-1)/sqrt(n[1]-2)),3)
#  
#  #n=10时
#  for (k in 1:m) {
#    x <- rexp(n[2],lammda)
#    lamstar<-numeric(B)
#    for(b in 1:B){
#      xstar <- sample(x,replace=TRUE)
#      lamstar[b]<-1/mean(xstar)
#    }
#    bias[k,2]<-mean(lamstar)-1/mean(x)
#    se[k,2]<-sd(lamstar)
#  }
#  round(c(bias2=mean(bias[,2]),se2=mean(se[,2]),bias2.theo=lammda/(n[2]-1),se2.theo=lammda*n[2]/(n[2]-1)/sqrt(n[2]-2)),3)
#  
#  #n=20时
#  for (k in 1:m) {
#    x <- rexp(n[3],lammda)
#    lamstar<-numeric(B)
#    for(b in 1:B){
#      xstar <- sample(x,replace=TRUE)
#      lamstar[b]<-1/mean(xstar)
#    }
#    bias[k,3]<-mean(lamstar)-1/mean(x)
#    se[k,3]<-sd(lamstar)
#  }
#  round(c(bias3=mean(bias[,3]),se3=mean(se[,3]),bias3.theo=lammda/(n[3]-1),se3.theo=lammda*n[3]/(n[3]-1)/sqrt(n[3]-2)),3)

## -----------------------------------------------------------------------------
## bias1        se1 bias1.theo   se1.theo 
## 0.564      2.074      0.500      1.443 
## bias2        se2 bias2.theo   se2.theo 
## 0.226      0.844      0.222      0.786 
## bias3        se3 bias3.theo   se3.theo 
## 0.103      0.500      0.105      0.496 

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12)
#  library(boot); library(MASS)
#  lammda<-2
#  n<-c(5,10,20);B<-1e3;m<-1e3
#  bias<-matrix(0, nrow =m, ncol = 3)
#  se<-matrix(0, nrow =m, ncol = 3)
#  
#  f.lam <- function(x,i) 1/mean(x[i])
#  
#  #n=5
#  for (k in 1:m) {
#    x <- rexp(n[1],lammda)
#    obj <- boot(data=x,statistic=f.lam,R=B)
#    bias[k,1]<-mean(obj$t)-obj$t0
#    se[k,1]<-sd(obj$t)
#  }
#  round(c(mean.bias=mean(bias[,1]),mean.se=mean(se[,1])),3)
#  
#  #n=10
#  for (k in 1:m) {
#    x <- rexp(n[2],lammda)
#    obj <- boot(data=x,statistic=f.lam,R=B)
#    bias[k,2]<-mean(obj$t)-obj$t0
#    se[k,2]<-sd(obj$t)
#  }
#  round(c(mean.bias=mean(bias[,2]),mean.se=mean(se[,2])),3)
#  
#  #n=20
#  for (k in 1:m) {
#    x <- rexp(n[3],lammda)
#    obj <- boot(data=x,statistic=f.lam,R=B)
#    bias[k,3]<-mean(obj$t)-obj$t0
#    se[k,3]<-sd(obj$t)
#  }
#  round(c(mean.bias=mean(bias[,3]),mean.se=mean(se[,3])),3)

## -----------------------------------------------------------------------------
## mean.bias   mean.se 
##     0.658     2.493
## mean.bias   mean.se 
##     0.224     0.835
## mean.bias   mean.se 
##     0.106     0.504

## -----------------------------------------------------------------------------
set.seed(111)
library(boot)
library(bootstrap) #for the law data

print(cor(law$LSAT, law$GPA))#直接统计的协方差

fff<-function(x){#统计量表达式
  return(cor(x[,1],x[,2]))
}


boot.t.ci <-#仿照书上的函数compute the bootstrap t CI
function(x, B = 500, R = 100, level = .95, statistic){
x <- as.matrix(x)#as.matrix()将x转换为矩阵形式
n <- nrow(x)
stat <- numeric(B); se <- numeric(B)

boot.se <- function(x, R, f) {#局部自举函数，计算对x自举的统计量f的标准差
x <- as.matrix(x); m <- nrow(x)
th <- replicate(R, expr = {
i <- sample(1:m, size = m, replace = TRUE)
f(x[i, ])
})#replicate重复执行R次语句expr
return(sd(th))
}

for (b in 1:B) {#正式开始，对x自举B次
j <- sample(1:n, size = n, replace = TRUE)
y <- x[j, ]#自举结果
stat[b] <- statistic(y)#theta.hat.(b)
se[b] <- boot.se(y, R = R, f = statistic)#使用局部自举函数计算theta.hat.(b)的标准差
}

stat0 <- statistic(x)#theta.hat
t.stats <- (stat - stat0) / se
se0 <- sd(stat)#theta.se
alpha <- 1 - level
Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)#寻找分位数
names(Qt) <- rev(names(Qt))#将名称进行反转,但不移动数值
CI <- rev(stat0 - Qt * se0)#置信区间
CI
}
boot.t.ci(law,statistic=fff)

## -----------------------------------------------------------------------------
set.seed(111)
library(bootstrap) #for the law data
x<-matrix(0,nrow=15,ncol=2)

print(cor(law$LSAT, law$GPA))#直接统计的协方差
x[,1]<-law$LSAT
x[,2]<-law$GPA

fff<-function(x,i){#统计量表达式
  return(cor(x[i,1],x[i,2]))
}
de<-boot(data=law,statistic=fff, R = 1e3)
boot.ci(de,type=c("all"))

## -----------------------------------------------------------------------------
library(boot)
set.seed(12)

x<-c(3,5,7,18,43,85,91,98,100,130,230,487)

boot.mean <- function(x,i) mean(x[i])
de <- boot(data=x,statistic=boot.mean, R = 1e3)
de
ci <- boot.ci(de,type=c("norm","basic","perc","bca"))
mean(x)

cat('norm =',ci$norm[2:3],'\nbasic =',ci$basic[4:5],'\nperc =',ci$percent[4:5],'\nBCa =',ci$bca[4:5])


## -----------------------------------------------------------------------------
n<-88
data(scor,package = 'bootstrap')
thetai<-numeric(n)

eigen<-eigen(cov(scor))$values
theta_hat<-eigen[5]/sum(eigen)

for (i in 1:n) {
  eigen<-eigen(cov(scor[-i,]))$values
  thetai[i]<-eigen[5]/sum(eigen)
}
cat('bias=',(n-1)*(mean(thetai)-theta_hat),"\n se=",
sqrt((n-1)*var(thetai)))


## -----------------------------------------------------------------------------
library('DAAG'); attach(ironslag)
set.seed(1)
n <- length(magnetic) #in DAAG ironslag
n
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1))
m<-1#指针
# leave-two-out cross validation
for (k in 2:n) {
  for(j in 1:(k-1)){#*很容易出错的k-1要带括号
    
    y <- magnetic[-c(k,j)]
    x <- chemical[-c(k,j)]
    
    #linear
    J1 <- lm(y ~ x)#拟合
    yhatk1 <- J1$coef[1] + J1$coef[2] * chemical[k]#预报k
    yhatj1 <- J1$coef[1] + J1$coef[2] * chemical[j]#预报j
    e1[m] <- magnetic[k] - yhatk1#预报误差k
    e1[m+1] <- magnetic[j] - yhatj1#预报误差j
    
    #quadratic
    J2 <- lm(y ~ x + I(x^2))
    yhatk2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
    J2$coef[3] * chemical[k]^2
    yhatj2 <- J2$coef[1] + J2$coef[2] * chemical[j] +
    J2$coef[3] * chemical[j]^2
    e2[m] <- magnetic[k] - yhatk2#预报误差k
    e2[m+1] <- magnetic[j] - yhatj2
    
    #exponential
    J3 <- lm(log(y) ~ x)
    logyhatk3 <- J3$coef[1] + J3$coef[2] * chemical[k]
    logyhatj3 <- J3$coef[1] + J3$coef[2] * chemical[j]
    yhatk3 <- exp(logyhatk3)
    yhatj3 <- exp(logyhatj3)
    e3[m] <- magnetic[k] - yhatk3#预报误差k
    e3[m+1] <- magnetic[j] - yhatj3
    
    #log-log
    J4 <- lm(log(y) ~ log(x))
    logyhatk4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
    logyhatj4 <- J4$coef[1] + J4$coef[2] * log(chemical[j])
    yhatk4 <- exp(logyhatk4)
    yhatj4 <- exp(logyhatj4)
    e4[m] <- magnetic[k] - yhatk4#预报误差k
    e4[m+1] <- magnetic[j] - yhatj4
    
    m<-m+2
  }
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))#平均平方预报误差
detach(ironslag)
#命令detach()结束使用数据集. 之前已经运行过attach().运行一次detach()只能删除上一次attach()的结果.所以通过多次运行detach()可以完全删除之前attach()的影响.再次运行attach()就不会报错了.


## -----------------------------------------------------------------------------
set.seed(1)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

m<-length(x);n<-length(y)

www<-function(data1,data2){##计算W2统计量的函数
  sum<-0
  f1<-ecdf(data1);f2<-ecdf(data2)
  m1<-length(data1);m2<-length(data2)
  for (i in 1:m1) {
    sum<-sum+(f1(data1[i])-f2(data1[i]))^2
  }
  for (i in 1:m2) {
    sum<-sum+(f1(data2[i])-f2(data2[i]))^2
  }
  sum<-sum*m1*m2/(m1+m2)^2
}

R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:(n+m)
reps <- numeric(R) #storage for replicates
t0 <- www(x,y)
for (i in 1:R) {
#generate indices k for the first sample
  k <- sample(K, size = 14, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  reps[i] <- www(x1,y1)
}
p <- mean(c(t0, reps) >= t0)
p

## -----------------------------------------------------------------------------
set.seed(12)
n1 <- 20;n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
m <- 1000;R<-999

maxout <- function(x, y) {#统计量函数maxout计算maximum number of extreme points
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
return(max(c(outx, outy)))
}

permu_count5<-function(x,y,R){#a permutation test based on maxout
  z<-c(x,y)
  resp<-numeric(R)
  t0<-maxout(x,y)
  K<-1:length(z)
  for (i in 1:R) {
  #generate indices k for the first sample
  k <- sample(K, size = length(x), replace = FALSE)
  x2 <- z[k]
  y2 <- z[-k] #complement of x1
  reps[i] <- maxout(x2,y2)
}
p <- mean(c(t0, reps) >= t0)
return(p)#返回p值(即t0所在分位数)
}

pvalue <- replicate(m, expr={#计算一型错误率
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
x <- x - mean(x) #centered by sample mean
y <- y - mean(y)
permu_count5(x, y,R)
})
mean(pvalue<0.05)


## -----------------------------------------------------------------------------
#设计函数
set.seed(1)
FFF<-function(N,b1,b2,b3,f00){
  x1<-rpois(N,1);x2<-rexp(N,1);x3<-rbinom(N,1,0.5)#随机量
  g <- function(alpha){#p为alpha函数
  p <- 1/(1+exp(-alpha-b1*x1-b2*x2-b3*x3))
  mean(p) - f00
  }
  #求根
  solution <- uniroot(g,c(-20,0))
  alpha <- solution$root
  return(alpha)
} 

## -----------------------------------------------------------------------------
out<-c(FFF(N=1e6,b1=0,b2=1,b3=-1,f00=0.1),
FFF(N=1e6,b1=0,b2=1,b3=-1,f00=0.01),
FFF(N=1e6,b1=0,b2=1,b3=-1,f00=0.001),
FFF(N=1e6,b1=0,b2=1,b3=-1,f00=0.0001))
round(out,3)

## -----------------------------------------------------------------------------
set.seed(1)
N<-1e2;b1<-0;b2<-1;b3<--1
x1<-rpois(N,1);x2<-rexp(N,1);x3<-rbinom(N,1,0.5)#随机量
s<-function(x){
  p <- 1/(1+exp(-x)*exp(-b1*x1-b2*x2-b3*x3))
  return(-log(mean(p)))
}
x<-seq(-20,10,length=200);y<-numeric(200)
for (i in 1:200) {
  y[i]<-s(x[i])
}
plot(x,y,type="l")#由于s函数中x可以为数组，所以直接curve(s)会报错

## -----------------------------------------------------------------------------
#生成这四个链
set.seed(111)
indept.Metropolis <- function(sigma_g, x0, N){#链生成
  #sigma_g为提议方差 x0为初始值 N链长 
    x = numeric(N)
    x[1] = x0
    u = runif(N)
    k = 0
    for(i in 2:N){
      y = rnorm(1, x[i-1], sigma_g)#产生候选值
      porb = exp(abs(x[i-1])-abs(y))#计算接受概率
      if(u[i] <= porb)
        x[i] = y  
      else{
        x[i] = x[i-1]
        k = k + 1
      }
    }
    return(list(x=x, k=k))
  }
chain_len = 5000
xt = numeric(chain_len)
sigma_g = c(0.05, 0.5, 2, 4)
x0 <- 3   # initial point
  indept1 <- indept.Metropolis(sigma_g[1], x0, chain_len)
  indept2 <- indept.Metropolis(sigma_g[2], x0, chain_len)
  indept3 <- indept.Metropolis(sigma_g[3], x0, chain_len)
  indept4 <- indept.Metropolis(sigma_g[4], x0, chain_len)

## -----------------------------------------------------------------------------
  refline <- qcauchy(c(.025, .975))
  indept <- cbind(indept1$x, indept2$x, indept3$x, indept4$x)
  for (j in 1:4) {
      plot(indept[,j], type="l",
           xlab=bquote(sigma == .(round(sigma_g[j],3))),
           ylab="X", ylim=range(indept[,j]))
      abline(h=refline)
  }

## -----------------------------------------------------------------------------
a <- ppoints(100)
QR <- -log(1-2*abs(a-1/2))*sign(a-1/2)#quantiles of Laplace
Q1 <- quantile(indept1$x[300:5000], a)
Q2 <- quantile(indept2$x[300:5000], a)
Q3 <- quantile(indept3$x[300:5000], a)
Q4 <- quantile(indept4$x[300:5000], a)
qqplot(QR, Q1, main="",xlab="Laplace Quantiles", ylab="Sample Quantiles")
abline(0,1,col='blue',lwd=2)
qqplot(QR, Q2, main="",xlab="Laplace Quantiles", ylab="Sample Quantiles")
abline(0,1,col='blue',lwd=2)
qqplot(QR, Q3, main="",xlab="Laplace Quantiles", ylab="Sample Quantiles")
abline(0,1,col='blue',lwd=2)
qqplot(QR, Q4, main="",xlab="Laplace Quantiles", ylab="Sample Quantiles")
abline(0,1,col='blue',lwd=2)

## -----------------------------------------------------------------------------
x<-seq(-20,20,length=200)
y<-1/2*exp(-abs(x))
hist(indept1$x,xlim=c(-20,20),breaks = 100,col = "pink",freq=FALSE)
lines(x,y,col="blue")
hist(indept2$x,xlim=c(-20,20),breaks = 100,col = "pink",freq=FALSE)
lines(x,y,col="blue")
hist(indept3$x,xlim=c(-20,20),breaks = 100,col = "pink",freq=FALSE)
lines(x,y,col="blue")
hist(indept4$x,xlim=c(-20,20),breaks = 100,col = "pink",freq=FALSE)
lines(x,y,col="blue")

## -----------------------------------------------------------------------------
print(c(indept1$k, indept2$k, indept3$k, indept4$k)/chain_len)#四个链的拒绝率

## ----echo=T-------------------------------------------------------------------
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
s1 <- sqrt(1-.9^2)
s2 <- sqrt(1-.9^2)
###### generate the chain #####
X[1, ] <- c(0, 0) #initialize
for (i in 2:N) {
x2 <- X[i-1, 2]
m1 <- .9 * x2
X[i, 1] <- rnorm(1, m1, s1)
x1 <- X[i, 1]
m2 <- .9* x1
X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]
x[1]


## ----echo=F-------------------------------------------------------------------
plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
lines(x[,2],col=2,lwd=2)
legend('bottomright',c(expression(X),expression(Y)),col=1:2,lwd=2)

## -----------------------------------------------------------------------------
lm(x[,2] ~ x[,1])
c(
cat('Means: ',round(colMeans(x),2)),
cat('\nStandard errors: ',round(apply(x,2,sd),2)),
cat('\nCorrelation coefficients: ', round(cor(x[,1],x[,2]),2)))

## -----------------------------------------------------------------------------
er<-x[,1]-0.899347*x[,1]+0.005348
qqnorm(er)
qqline(er, col = "red", lwd = 2, lty = 2)

hist(er,xlim=c(-.5,.5),breaks = 30,col = "pink",freq=FALSE)
xxx<-seq(-.5,.5,length=200);yyy<-dnorm(xxx,mean(er),sd(er))
lines(xxx,yyy,col="blue",lwd=2)

## -----------------------------------------------------------------------------
ks.test(scale(er),"pnorm")#scale(er)要先对数据进行标准化

## -----------------------------------------------------------------------------
wit<-100
num<-1:round(length(er)/wit)
sd<-numeric(length(num))
for (i in num) {
  sd[i]<-sd(er[(wit*(i-1)+1):(wit*i)])
}
plot(num,sd)

## -----------------------------------------------------------------------------
rm(list = ls())#删除当前工作环境中的所有对象。具体而言，这会删除所有已分配的变量、数据框、函数等，使得工作环境中不再包含这些对象

set.seed(123)
library("coda")
Gelman.Rubin <- function(psi) {#计算GR统计量
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
  psi <- as.matrix(psi)#转为矩阵
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic!!!!!!!!
  return(r.hat)
}
f <- function(x, sigma) {
  if (any(x < 0)) return (0)
  stopifnot(sigma > 0)
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}
gen.chain<-function(sigma,N,x1){
  x<-numeric(N)
  x[1]<-x1
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den) x[i] <- y 
    else x[i] <- xt
  }
  return(x)
}
n <- 10000
sigma <- 4
k<-3  #number of chains
b <- 1000 #burn-in length

#choose overdispersed initial values
x0 <- c(rchisq(1, df=1),rchisq(1, df=2),rchisq(1, df=3))
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k)
  X[i, ] <- gen.chain(sigma, n, x0[i])

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

## -----------------------------------------------------------------------------
print(Gelman.Rubin(psi))

## -----------------------------------------------------------------------------
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l",xlab=i, ylab=bquote(psi))

## -----------------------------------------------------------------------------
#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)
abline(h=1.2,lty=1)

## -----------------------------------------------------------------------------
l1<-as.mcmc(X[1,])
l2<-as.mcmc(X[2,])
l3<-as.mcmc(X[3,])
l<-mcmc.list(l1,l2,l3)
gelman.diag(l)

## -----------------------------------------------------------------------------
gelman.plot(l)

## -----------------------------------------------------------------------------
#区间
u<-c(11,8,27,13,16,0,23,10,24,2)
v<-c(12,9,28,14,17,1,24,11,25,3)
#对数似然的导数
gl<-function(lambda){
  sum((v*exp(-lambda*v)-u*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v)))
}
#直接求解MLE
lam_mle<-uniroot(gl,c(0,1),extendInt = "yes")$root

#EM算法:
lam0<-0.05#初始化
N<-1e3
tol<-1e-6
for(iter in 1:N){
  lam1<-1/(1/lam0-1/length(u)*gl(lam0))
  if((abs(lam1-lam0)/lam0)<=tol) break
  lam0<-lam1
}
lam_EM<-lam1

#结果
c("MLE"=lam_mle,"EM"=lam_EM)

## -----------------------------------------------------------------------------
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}

## -----------------------------------------------------------------------------
#enter the payoff matrix
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
s <- solve.game(A+2)

#结果
round(cbind(s$x, s$y), 7)

## -----------------------------------------------------------------------------
my_list <- list(1, "a", TRUE)
a<-unlist(my_list)
a;is.atomic(a)

## -----------------------------------------------------------------------------
my_list <- list(1, "a", TRUE)
a<-as.vector(my_list)
str(a);is.list(a)

## -----------------------------------------------------------------------------
my_vector <- c(1, 2, 3, 4, 5,list(1,"1"))
dim(my_vector)

## -----------------------------------------------------------------------------
# 创建一个数据框
my_df <- data.frame(
  A = c(1, 2, 3),
  B = c("a", "b", "c"),
  C = c(TRUE, FALSE, TRUE)
)
# 使用as.matrix()将数据框转换为矩阵
as.matrix(my_df)

## -----------------------------------------------------------------------------
em30 <- data.frame(matrix(nrow=3, ncol=0),row.names = c("p","q","t"))
dim(em30);em30

em03 <- data.frame(matrix(nrow=0, ncol=3))#无col.names =,否则会被当成一列数据^
colnames(em03) <- c("Column1", "Column2", "Column3")
dim(em03);em03

## -----------------------------------------------------------------------------
a <- data.frame()
dim(a);a

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)#计算向量 x的最小值和最大值。na.rm = TRUE 表示在计算范围时要移除缺失值（NA）
(x - rng[1]) / (rng[2] - rng[1])
}

scale01(c(1,2,3,4,5))

## -----------------------------------------------------------------------------
scale011 <- function(x) {
if (is.numeric(x)==F) return(x)
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
# 创建一个数据框
my_df <- data.frame(
  A = c(1, 2, 3,NA),
  B = c("a", "b", "c","p"),
  C = c(4, 5, 6,7),
  D = c(7, 8, 9,9)
)
data.frame(lapply(my_df, scale011))

## ----eval=FALSE---------------------------------------------------------------
#  library(dplyr)
#  # 使用mutate_if()用于对数据框的某些列进行变换。条件为is.numeric为真
#  mutate_if(my_df,is.numeric,scale01)#scale01(x)会报错

## -----------------------------------------------------------------------------
num<-sapply(my_df,is.numeric)
b<-which(num>0)#选择数值型的列索引
for (i in b)
  my_df[,i]<-scale01(my_df[,i])
my_df

## -----------------------------------------------------------------------------
df <- data.frame(
  A = c(1, 2, 3,NA),
  C = c(4, 5, 6,7),
  D = c(7, 8, 9,9))
vapply(df,sd,FUN.VALUE=c(sd=0),na.rm = TRUE)

## -----------------------------------------------------------------------------
my_df <- data.frame(
  A = c(1, 2, 3,NA),
  B = c("a", "b", "c","p"),
  C = c(4, 5, 6,7),
  D = c(7, 8, 9,9))

numeric_columns <- names(my_df)[sapply(my_df, is.numeric)]

vapply(my_df[numeric_columns],sd, numeric(1),na.rm = TRUE)

## -----------------------------------------------------------------------------
    library(Rcpp)
    dir_cpp <- './' #当前文件夹
    source(paste0(dir_cpp,'gibbsR.R')) #paste0(dir_cpp, "meanC.cpp") 的作用是将两个字符向量连接起来，生成一个新的字符向量。
    sourceCpp(paste0(dir_cpp,'gibbsC.cpp'))
    #sourceCpp()就是加载括号内地址cpp文件
    library(microbenchmark)
    n=5;a=1;b=1
    ts <- microbenchmark(gibbR=gibbsR(100,10,a,b,n), gibbC=gibbsC(100,10,a,b,n))
    
    knitr::kable(summary(ts)[,c(1,3,5,6)],format='markdown')

## -----------------------------------------------------------------------------
c(3290.55	,4055.40,	4616.55	)/c(408.20,	525.15,	561.90	)

