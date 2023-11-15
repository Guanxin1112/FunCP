#' Simulated data from the Subgroup learning in functional regression under the RKHS framework
#'
#' The function for simulating data from the Subgroup learning in functional regression under the RKHS framework
#'
#' @rdname simdata
#' @param n sample size.
#' @param m functional individual size.
#' @param oo the subset vector
#' @param kenel_sigma the kernel bandwidth of a gaussian kernel.
#' @param gamma the threshold parameter for default.
#'
#' @return A list with the elements
#' \item{y}{The response variable.}
#' \item{x}{The scalar covariate with threshold.}
#' \item{z}{A vector of covariates.}
#' \item{s}{A vector of locations or times.}
#' \item{k}{A matrix of kernel gram.}
#'
#' @author Xin Guan, Xu Liu and Jinhong You
#' @keywords simdata

#' @importFrom MASS mvrnorm
#' @importFrom stats runif rnorm
#' @export
#'
#' @examples
#'
#' ## simulated data
#' ptm <- proc.time()
#' n <- 100
#' m <- 10
#' kenel_sigma <- 0.2
#' oo <- c(2, 3)
#' gamma <- c(-1,1)
#' dat <- simdata(n,m,kenel_sigma,oo,gamma)
#' head(dat)
#' proc.time() - ptm

### generate the data
simdata<-function(n,m,kenel_sigma,oo,gamma){

  beta1f<-function(s) (1-s)^3
  beta2f<-function(s) exp(-s^2)
  #beta3f<-function(s) sin(2*pi*s) + s^2
  beta3f<-function(s) sin(pi*s) + s^3

  delta1f<-function(s) (1-s)^2
  #delta2f<-function(s) 2*s*(1-s) + cos(2*pi*s)
  delta2f<-function(s) exp(-5*s)


  #individual function
  phi1<-function(s) sqrt(2)*sin(2*pi*s)
  phi2<-function(s) sqrt(2)*cos(2*pi*s)

  ### Reproducing Kernel function
  gaussian<-function(x,kenel_sigma){
    exp(-x^2/(2*kenel_sigma^2))
  }


  #### Function to obtain the Gram Matrix
  getgram<-function(s,kenel_sigma){
    m = length(s)
    Gram = matrix(0,m,m)
    for(i in 1:m){
      for(j in 1:m){
        Gram[i,j] = gaussian(s[i]-s[j],kenel_sigma)
      }
    }
    return(Gram)
  }


  Mean = c(0,0)
  #Sigma = matrix(c(1,0.5,0.5^2,0.5,1,0.5,0.5^2,0.5,1),3,3)
  Sigma = matrix(c(1,0.5,0.5,1),2,2)
  X1 = mvrnorm(n, Mean, Sigma)
  X = cbind(rep(1,n), X1)

  S = runif(m,0,1)
  beta = cbind(beta1f(S),beta2f(S),beta3f(S))
  delta = cbind(delta1f(S),delta2f(S))

  K = getgram(S, kenel_sigma) # kernel gaussian gram

  ###change variable
  tX = X[,oo]
  Z1 = rnorm(n, 0, 1)
  Z2 = rnorm(n, 1, 1)#runif(n,0,1)
  Z = cbind(Z1, rep(1,n), Z2)

  ind = (Z[,1] + Z[,-1]%*%gamma >0)

  Y = matrix(0,n,m)
  eta = matrix(0,n,m)
  for(i in 1:n){
    eps = rnorm(m,0,sqrt(0.1))
    eta[i,] = rnorm(1,0,1)*phi1(S) + rnorm(1,0,0.5^2)*phi2(S)
    tem1 = X[i,1]*beta[,1] + X[i,2]*beta[,2]+ X[i,3]*beta[,3]
    tem2 = tX[i,1]*delta[,1] + tX[i,2]*delta[,2]
    Y[i,] =tem1 + tem2 * as.numeric(ind[i]) + eta[i,] + eps
  }

  return(list(X=X,Y=Y,S=S,Z=Z,K=K))
}





