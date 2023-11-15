#' Simulated data from the Subgroup learning in functional regression under the RKHS framework
#'
#' The function for simulating data from the Subgroup learning in functional regression under the RKHS framework
#'
#' @rdname smooth_ls
#' @param Y A vector of response.
#' @param X A vector of covariates
#' @param Z A scalar covariate with threshold.
#' @param K A matrix of kernel gram.
#' @param oo the subset vector
#' @param lambda the parameter.
#' @param gamma.ini the threshold parameter for initial.
#' @param h the smooth parameter
#'
#' @return A list with the elements
#' \item{res}{The residuals.}
#' \item{b}{The regression patameters}
#' \item{Yhat}{The estimated response}
#' \item{etahat}{The estimated eta function.}
#' \item{epshat}{The estimated eps function.}
#'
#' @author Xin Guan, Xu Liu and Jinhong You
#' @keywords smooth_ls

#' @importFrom MASS mvrnorm
#' @importFrom stats runif rnorm dnorm pnorm optim
#' @importFrom Matrix bdiag
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
#' data <- simdata(n,m,kenel_sigma,oo,gamma)
#' X = data$X; Y=data$Y; S=data$S; K=data$K;Z=data$Z
#' gamma.ini = c(0,0)
#' lambda = 0.01
#' h = log(n)*n^(-0.5)
#' result <- smooth_ls(X, Y, Z, K, lambda, oo, gamma.ini, h)
#' proc.time() - ptm

### estimation without considering covariance
smooth_ls<-function(X, Y, Z, K,lambda,oo,gamma.ini,h){
  n = dim(X)[1]
  m = dim(K)[1]
  iter = 0
  tol = 10


  Gsmooth<-function(xx, h){
    #kerf<-function(s)  {(1/sqrt(2*pi)) * exp(-s^2/2)}
    #integrate(kerf, lower =-1000, upper = xx/h, abs.tol = 0L)$value
    pnorm(xx/h) + (xx/h)*dnorm(xx/h)
  }



  penest_full<-function(X, Y, Z, K, lambda,oo,gamma,h){
    n = dim(X)[1]
    p = dim(X)[2]
    m = dim(K)[1]
    d = length(oo)

    ind = Z[,1] + Z[,-1] %*% gamma

    G = apply(ind,1,function(xx)return(Gsmooth(xx, h)))

    omega = as.matrix(bdiag(diag(1,p,p),diag(1,d,d)))


    D1 = array(0,dim=c(m*(p+d),m*(p+d),n))
    D2 = array(0,dim=c(m*(p+d),1,n))
    for(i in 1:n){
      temp = c(X[i,],X[i,oo]*G[i])
      N = kronecker(t(temp),K,FUN = "*")
      D1[,,i] = t(N) %*% N #+ m*lambda * kronecker(omega,K,FUN = "*")
      D2[,,i] = t(N) %*% Y[i,]
    }
    #b = ginv(apply(D1,c(1,2),sum)) %*% apply(D2,c(1,2),sum)
    D1_temp = apply(D1,c(1,2),sum) + n*m*lambda*kronecker(omega,K,FUN = "*")
    b = solve(D1_temp,tol=1e-20) %*% apply(D2,c(1,2),sum)
    #system is computational singular reciprocal condition number = 2e-17

    beta1fhat = K %*% b[1:m]
    beta2fhat = K %*% b[(m+1):(2*m)]
    beta3fhat = K %*% b[(2*m+1):(3*m)]

    delta1fhat = K %*% b[(3*m+1):(4*m)]
    delta2fhat = K %*% b[(4*m+1):(5*m)]

    Yhat = res = etahat = epshat = matrix(0,n,m)
    for(i in 1:n){
      temp = c(X[i,],X[i,oo]*G[i])
      N = kronecker(t(temp),K,FUN = "*")
      Yhat[i,] = N %*% b
      res[i,] = Y[i,]-Yhat[i,]
      d = solve(t(K) %*% K,tol=1e-20) %*% t(K) %*% as.vector(res[i,])
      etahat[i,] = K %*% d
      epshat[i,] = res[i,] - etahat[i,]
    }

    # res = Y-Yhat
    return(list(b=b,res=res,Yhat=Yhat,etahat=etahat,
                epshat=epshat,beta1fhat=beta1fhat,
                beta2fhat=beta2fhat,beta3fhat=beta3fhat,
                delta1fhat=delta1fhat,delta2fhat=delta2fhat))
  }



  while(iter<20 & tol>1e-6){
    resultini = penest_full(X, Y, Z, K, lambda,oo,gamma=gamma.ini,h)
    b = resultini$b

    fn<-function(xx){

      ind =  Z[,1] + Z[,-1] %*% xx

      G = apply(ind,1,function(s)return(Gsmooth(s, h)))

      Yhat = matrix(0,n,m)
      for(i in 1:n){
        temp = c(X[i,],X[i,oo]*G[i])
        N = kronecker(t(temp),K,FUN = "*")
        Yhat[i,] = N %*% b
        # Yhat[i,] = Y[i,]-betaest %*% X[i,]-(deltaest %*% X[i,oo]) * G[i]
      }

      yvec = as.vector(t(Y))
      yhatvec = as.vector(t(Yhat))
      rss = t(yvec-yhatvec) %*% (yvec-yhatvec)

      return(rss)
    }

    R = optim(par=gamma.ini, fn=fn,gr=NULL)
    gammanew = R$par

    tol = max(abs(gammanew-gamma.ini))
    iter = iter + 1

    gamma.ini = gammanew

  }

  #resultini = penest_full(X, Y, Z, K, lambda,oo,gamma=gammanew,h)

  return(list(res=resultini$res,b=resultini$b,Yhat = resultini$Yhat,
              etahat=resultini$etahat,epshat=resultini$epshat,
              gammanew=gammanew,beta1fhat=resultini$beta1fhat,
              beta2fhat=resultini$beta2fhat,beta3fhat=resultini$beta3fhat,
              delta1fhat=resultini$delta1fhat,delta2fhat=resultini$delta2fhat))

}
