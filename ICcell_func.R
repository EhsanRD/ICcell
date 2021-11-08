#rm(list = ls())

# Functions for generating design matrices
SWdesmat <- function(T) {
  Xsw <- matrix(data=0, ncol = T, nrow = (T-1))
  for(i in 1:(T-1)) {
    Xsw[i,(i+1):T] <- 1
  }
  return(Xsw)
}

#Designs missing cluster-period cells:

###this calculates the variance of the treatment effect estimator for the given inputs:###

#General function (taken from ExpDecayVar.r in Information content folder)
CRTVarGeneral <- function(Xmat, m, rho0, r, type) {

  totalvar <- 1
  sig2CP <- rho0*totalvar
  sig2E <- totalvar - sig2CP

  sig2 <- sig2E/m

  T <- ncol(Xmat)
  K <- nrow(Xmat)

  Xvec <- as.vector(t(Xmat))
  stackI <- matrix(rep(diag(1,T)), nrow=K*T, ncol=T, byrow=TRUE)
  Zmat <- cbind(stackI[!is.na(Xvec),], Xvec[!is.na(Xvec)])

  #Variance matrix for one cluster, with decay in correlation over time
  #Vi <- diag(sig2,T) + matrix(sig2CP,nrow=T, ncol=T)
  #Constant decay var if type==0
  if(type==0) { 
    Vi <-diag(sig2 +(1-r)*sig2CP, T) + matrix(data=sig2CP*r, nrow=T, ncol=T)
  }
  #exponential decay structure
  if(type==1) { 
    Vi <- diag(sig2,T) + sig2CP*(r^abs(matrix(1:T,nrow=T, ncol=T, byrow=FALSE) - matrix(1:T,nrow=T, ncol=T, byrow=TRUE)))
  }
  #Variance matrix for all clusters
  Vall <- kronecker(diag(1,K), Vi)
  Vall <- Vall[!is.na(Xvec),!is.na(Xvec)]

  #vartheta <- solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)]

  #xtry <- try(solve((t(Zmat)%*%solve(Vall)%*%Zmat))) 
  #if('try-error' %in% class(xtry)) return(NA)

  #there will be problems if Zmat is not of full column rank
  if(rankMatrix(Zmat)[1] < ncol(Zmat)) return(NA)
  else return(solve((t(Zmat)%*%solve(Vall)%*%Zmat))[ncol(Zmat),ncol(Zmat)])

  #return(vartheta)
}

