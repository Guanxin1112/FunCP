#' Simulated data from the Subgroup learning in functional regression under the RKHS framework
#'
#' The function for simulating data from the Subgroup learning in functional regression under the RKHS framework
#'
#' @rdname smooth_weightls
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
#' @keywords smooth_weightls

#' @importFrom MASS mvrnorm ginv
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
#' result <- smooth_weightls(X, Y, Z, K, lambda, oo, gamma.ini, h)
#' proc.time() - ptm

### estimation by considering covariance
smooth_weightls<-function(X,Y, Z, K,lambda,oo,gamma.ini,h){

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


  #### Function to obtain weighted estimation
  weight_penest_full<-function(X,Y, Z,K,varPhi.inv,lambda,oo,gamma,h){
    n = dim(X)[1]
    p = dim(X)[2]
    m = dim(K)[1]
    d = length(oo)

    ind = Z[,1] + Z[,-1] %*% gamma

    G = apply(ind,1,function(xx)return(Gsmooth(xx, h)))
    # plot(sort(as.numeric((ind>0))))
    # lines(sort(G))

    omega = as.matrix(bdiag(diag(1,p,p),diag(1,d,d)))

    D1 = array(0,dim=c(m*(p+d),m*(p+d),n))
    D2 = array(0,dim=c(m*(p+d),1,n))
    for(i in 1:n){
      temp = c(X[i,],X[i,oo]*G[i])
      #temp = c(X[i,],X[i,oo]*as.numeric(ind[i]>0))
      N = kronecker(t(temp),K,FUN = "*")
      D1[,,i] = t(N) %*% varPhi.inv %*% N #+ m*lambda * kronecker(omega,K,FUN = "*")
      D2[,,i] = t(N) %*% varPhi.inv %*% Y[i,]
    }
    #weight_b = ginv(apply(D1,c(1,2),sum)) %*% apply(D2,c(1,2),sum)
    D1_temp = apply(D1,c(1,2),sum) + n*m*lambda*kronecker(omega,K,FUN = "*")
    weight_b = solve(apply(D1_temp,c(1,2),sum),tol=1e-20) %*% apply(D2,c(1,2),sum)


    beta1fhat = K %*% weight_b[1:m]
    beta2fhat = K %*% weight_b[(m+1):(2*m)]
    beta3fhat = K %*% weight_b[(2*m+1):(3*m)]

    delta1fhat = K %*% weight_b[(3*m+1):(4*m)]
    delta2fhat = K %*% weight_b[(4*m+1):(5*m)]


    Yhat = res = etahat = epshat = matrix(0,n,m)
    rss = 0
    tr_A = 0
    for(i in 1:n){
      temp = c(X[i,],X[i,oo]*G[i])
      N = kronecker(t(temp),K,FUN = "*")
      Yhat[i,] = N %*% weight_b
      res[i,] = Y[i,]-Yhat[i,]
      d = solve(t(K) %*% K,tol=1e-20) %*% t(K) %*% as.vector(res[i,])
      etahat[i,] = K %*% d
      epshat[i,] = res[i,] - etahat[i,]
      # rss = rss + sum((Y[i,] - Yhat[i,])^2)
      # A = Z %*% ginv(apply(Z1,c(1,2),sum)) %*% t(Z) %*% ginv(varPhi)
      # tr_A = tr_A + sum(diag(A))
    }
    # gcv = (eps*(m*n)^(-1))/((1-tr_A*(n*m)^(-1))^2)



    # yvec = as.vector(t(Y))
    # yhatvec = as.vector(t(Yhat))
    # rss = t(yvec-yhatvec) %*% (yvec-yhatvec)
    #
    # Zmat = kronecker(X,K,FUN = "*")
    # varPhimat = kronecker(diag(1,n,n),varPhi,FUN = "*")
    # Mmat = ginv(varPhimat)
    # temp = ginv(t(Zmat) %*% Mmat %*% Zmat + n*m*lambda * kronecker(diag(1,p,p),K,FUN = "*"))
    # A = Zmat %*% temp %*% t(Zmat) %*% Mmat
    # tr_A1 = sum(diag(A))
    # gcv = rss*(n*m)^(-1)/((1-tr_A*(n*m)^(-1))^2)


    return(list(Yhat=Yhat,res=res,etahat=etahat,
                epshat=epshat, weight_b=weight_b,
                beta1fhat=beta1fhat,beta2fhat=beta2fhat,beta3fhat=beta3fhat,
                delta1fhat=delta1fhat,delta2fhat=delta2fhat))
  }


  ###Function to obtain variance
  varest<-function(res, K){
    n = dim(res)[1]
    m = dim(K)[1]
    #resmat = matrix(res,n,m,byrow = T) # ATTENTION !
    resmat = matrix(res,n,m,byrow = F)
    etahat = matrix(0,n,m)
    epshat = matrix(0,n,m)
    Rhat = matrix(0,m,m)
    for(i in 1:n){
      d = ginv(t(K) %*% K) %*% t(K) %*% as.vector(resmat[i,])
      etahat[i,] = K %*% d
      Rhat = Rhat + (n^(-1)) * etahat[i,] %*% t(etahat[i,])
      epshat[i,] = (resmat[i,] - etahat[i,])^2
    }
    epstemp = apply(epshat,2,sum)
    a = ginv(t(K) %*% K,tol=1e-20) %*% t(K) %*% (epstemp*n^(-1))
    epssigma = K %*% a
    epssigma = mean(as.vector(epssigma)) # ATTENTION !

    # spectral decomposition
    eigenR = eigen(Rhat)
    numR  = sum(eigenR$values > 1e-15)
    eigenRval = eigenR$values[1L:numR]
    eigenRvec = eigenR$vectors[,1L:numR,drop=FALSE]
    # V = eigenRvec %*% diag(sqrt(eigenRval))
    etacov = eigenRvec %*% diag(eigenRval) %*% t(eigenRvec)

    #varPhi = etacov + diag(as.vector(epssigma,m,m)
    varPhi = etacov + diag(epssigma,m)

    return(varPhi)
  }



  n = dim(X)[1]
  m = dim(K)[1]
  iter = 0
  toltem = 10
  while(iter<20 & toltem>1e-6){
    resultini = penest_full(X, Y, Z, K, lambda,oo,gamma=gamma.ini,h)
    res = resultini$res
    varPhi = varest(res=res, K)
    varPhi.inv = solve(varPhi,tol=1e-20)

    result = weight_penest_full(X,Y,Z,K,varPhi.inv,lambda,oo,gamma=gamma.ini,h)
    weight_b = result$weight_b

    fn<-function(xx){

      ind =  Z[,1] + Z[,-1] %*% xx

      G = apply(ind,1,function(s)return(Gsmooth(s, h)))

      Yhat = matrix(0,n,m)
      for(i in 1:n){
        temp = c(X[i,],X[i,oo]*G[i])
        N = kronecker(t(temp),K,FUN = "*")
        Yhat[i,] = N %*% weight_b
        # Yhat[i,] = Y[i,]-betaest %*% X[i,]-(deltaest %*% X[i,oo]) * G[i]
      }

      tem = Y - Yhat
      rss = 0
      for(i in 1:n){rss = rss +  t(tem[i,]) %*% varPhi.inv %*% tem[i,]}

      # yvec = as.vector(t(Y))
      # yhatvec = as.vector(t(Yhat))
      # rss1 = t(yvec-yhatvec) %*% (yvec-yhatvec)

      return(rss)
    }

    R = optim(par=gamma.ini, fn=fn,gr=NULL)
    gammanew = R$par

    toltem = max(abs(gammanew-gamma.ini))
    iter = iter + 1

    gamma.ini = gammanew

  }

  # resultini = penest_full(X, Y, Z, K, lambda,oo,gamma=gammanew,h)
  # res = resultini$res
  # varPhi = varest(res=res, K)
  # varPhi.inv = ginv(varPhi)
  #
  # result = weight_penest_full(X,Y,Z,K,varPhi.inv,lambda,oo,gamma=gammanew,h)
  # weight_b = result$weight_b


  return(list(b=weight_b,res=result$res,Yhat=result$Yhat,
              etahat=result$etahat,epshat=result$epshat,
              gammanew=gammanew,beta1fhat=result$beta1fhat,
              beta2fhat=result$beta2fhat,beta3fhat=result$beta3fhat,
              delta1fhat=result$delta1fhat,delta2fhat=result$delta2fhat))

}
