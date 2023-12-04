gibbsR <- function(N, thin,a,b,n) {
  mat <- matrix(nrow = N, ncol = 2)#x,y分别抽N次的结果
  x <- y <- 0#初始值为0
  for (i in 1:N) {
    for (j in 1:thin) {#不直接抽样而是迭代几次 相当于间隔取样,降低结果相关性
      x <- rbinom(1,n,y)
      y <- rbeta(1, x + a, n-x+b)
    }
    mat[i, ] <- c(x, y)
  }
  mat
}
