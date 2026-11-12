hazard2D.LL<-function(t.grid,z.grid,O.tz,E.tz,bs.grid,cv=TRUE)

{
  # Univariate kernel: sextic kernel
  K<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}

  ## 1. Compute kernel evaluations and moments of the kernel
  ## 1.1. Kernel evaluations at dimension t
  M<-length(t.grid)
  Tt<-matrix(rep(t.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
  # Tt is an MxM matrix with the values of the grid
  Tx<-Tt-t(Tt)
  nb<-nrow(bs.grid)
  K.bt<-array(0,dim=c(M,M,nb))
  # Each element of K.bt is the Kernel evaluated at ut below
  ut<-array(0,dim=c(M,M,nb))
  for(b in 1:nb)
  {
    ut[,,b]=Tx/bs.grid[b,1]
    K.bt[,,b]<-apply(ut[,,b],1:2,K)
    # This works and returns the value of the Kernel function at each individual value of the array ut
    K.bt[,,b]<-K.bt[,,b]/(bs.grid[b,1])
  }

  ## 1.2. Dimension z (evaluation points are the same but bandwidth not necessarily)
  N <- length(z.grid)
  Zt <- matrix(rep(z.grid, times=N),nrow = N, ncol = N,byrow=FALSE)
  # Zt is an NxN matrix with the values of the grid
  Zx<-Zt-t(Zt)
  K.bz<-array(0,dim=c(N,N,nb))
  # Each element of K.bz is the Kernel evaluated at uz below
  uz<-array(0,dim=c(N,N,nb))
  for(b in 1:nb)
  {
    uz[,,b]=Zx/bs.grid[b,2]
    K.bz[,,b]<-apply(uz[,,b],1:2,K)
    # This works and returns the value of the Kernel function at each individual value of the array ut
    K.bz[,,b]<-K.bz[,,b]/(bs.grid[b,2])
  }
  rm(ut,uz,Tt,Zt)

  ## 1.3. Moment calculations involved in the estimator
  c0<-c10<-c11<-d00<-d01<-d11<-det.D<-den<-array(0,dim=c(M,N,nb))
  for (b in 1:nb)
  {
    c0[,,b]<-(K.bt[,,b]) %*% (E.tz) %*%t( K.bz[,,b])
    c10[,,b]<-(K.bt[,,b] * Tx) %*% (E.tz) %*%t(K.bz[,,b])
    c11[,,b]<-(K.bt[,,b]) %*% (E.tz) %*% t(K.bz[,,b]*Zx)
    d00[,,b]<-(K.bt[,,b] * (Tx^2)) %*% (E.tz) %*% t(K.bz[,,b])
    d01[,,b]<-(K.bt[,,b] * Tx) %*% (E.tz) %*% t(K.bz[,,b]*Zx)
    d11[,,b]<-(K.bt[,,b]) %*% (E.tz) %*% t(K.bz[,,b]*(Zx^2))
    det.D[,,b]<- d00[,,b]*d11[,,b]-d01[,,b]^2
    den[,,b]<-det.D[,,b]*c0[,,b] -
      ( c10[,,b]*(d11[,,b]*c10[,,b]-d01[,,b]*c11[,,b])
        + c11[,,b]*(d00[,,b]*c11[,,b] - d01[,,b]*c10[,,b]) )
  }

  ## 2. Create a function to compute the LL hazard from results above
  ##    and for one single bandwidth
  hazard.b<-function(K.bt,K.bz,Tx,Zx,b,O.tz,c10,c11,d00,d01,d11,det.D,den)
  {
    num1<-((K.bt[,,b]) %*% (O.tz)  %*%  t(K.bz[,,b])) *(det.D[,,b])
    num2<-( (K.bt[,,b] * Tx) %*% (O.tz)  %*%  t(K.bz[,,b]))*(d11[,,b]*c10[,,b] - d01[,,b]*c11[,,b])
    num3<-( (K.bt[,,b]) %*% (O.tz)  %*%  t(K.bz[,,b]*Zx))*(d00[,,b]*c11[,,b] - d01[,,b]*c10[,,b])
    num<-num1-num2-num3
    estim<-num/den[,,b] ### alpha(t,z)
    estim[den[,,b]==0]<-NA
    estim[estim<0]<-NA
    return(estim)
  }

  cv.values <- NULL

  ## 3. Compute now the estimator with possible CV-bandwidth
  if(cv==TRUE & nb>1)
  {
    ## 3.1. The leave-one-out (loo) LL-hazard estimator
    # Create a function similar to hazard.b but leave-one-out
    hazard.loo.b<-function(K.bt,K.bz,Tx,Zx,b,O.tz,c10,c11,d00,d01,d11,
                           det.D,den,M,N)
    {
      estim.ij<-double(M*N)
      ij<-l<-0
      for (i in 1:M)
      {
        for (j in 1:N)
        {
          l<-l+1
          if (O.tz[i,j]>0) {
            Oij.prev<-O.tz[i,j];O.tz[i,j]<- O.tz[i,j]-1;ij=1
          }
          num1<-( (K.bt[i,,b]) %*% (O.tz)  %*%  (K.bz[j,,b])) *(det.D[i,j,b])
          num2<-( (K.bt[i,,b] * Tx[i,]) %*% (O.tz)  %*%  K.bz[j,,b])*(d11[i,j,b]*c10[i,j,b] - d01[i,j,b]*c11[i,j,b])
          num3<-( (K.bt[i,,b]) %*% (O.tz)  %*%  (K.bz[j,,b]*Zx[j,]) )*(d00[i,j,b]*c11[i,j,b] - d01[i,j,b]*c10[i,j,b])
          num<-num1-num2-num3
          estim<-num/den[i,j,b]; estim[den[i,j,b]==0]<-NA
          estim.ij[l]<-estim
          if (ij>0) {O.tz[i,j]<-Oij.prev; ij=0}
        }
      }
      return(estim.ij)
    }

    # Finally evaluate the above function above for each bandwidth in the grid
    estim.loo.tz.bs <- sapply(1:nb, function(b){
      loo_vec <- hazard.loo.b(K.bt,K.bz,Tx,Zx,b,O.tz,c10,c11,d00,d01,d11,det.D,den,M,N)
      loo_mat <- matrix(loo_vec, nrow = M, ncol = N, byrow = TRUE)
      return(as.double(loo_mat))
    })
    # estim.loo.tz.bs is a matrix with Mt*Mz rows and nb columns

    ## 3.2. The usual estimator (not loo)
    estim.tz.bs<-sapply(1:nb, function(b){
      return(hazard.b(K.bt,K.bz,Tx,Zx,b,O.tz,c10,
                      c11,d00,d01,d11,det.D,den))})
    ## estim.tz.bs is a matrix with Mt*Mz rows and nb columns

    ## 3.3. Compute the CV-banwidth choice with the above matrices
    vec.E.tz<-as.double(E.tz)
    vec.O.tz<-as.double(O.tz)
    #  The cross-validation function
    CV.b<-function(b)
    {
      dif1.b<-sum((estim.tz.bs[,b])^2 * vec.E.tz,na.rm=TRUE)
      dif2.b<-sum(estim.loo.tz.bs[,b]* (vec.O.tz),na.rm=TRUE)
      cv.b<-dif1.b-2*dif2.b
      if (cv.b==Inf | cv.b==0) cv.b<-NA
      return(cv.b)
    }

    cv.values <- sapply(1:nb, CV.b)
    bcv<-which.min(cv.values)
    hi.tz<-estim.tz.bs[,bcv]
    hi.tz<-matrix(hi.tz,M,N,byrow=FALSE)
    hi.tz[which(hi.tz<0)]<-NA

  } else {

    estim.tz.b<-hazard.b(K.bt,K.bz,Tx,Zx,b=1,O.tz,c10,
                         c11,d00,d01,d11,det.D,den)
    hi.tz<-matrix(as.double(estim.tz.b),M,N,byrow=F)
    hi.tz[which(hi.tz<0)]<-NA

    bcv <- 1

  }

  return(list(hi.tz=hi.tz, bcv=bs.grid[bcv,], cv.values=cv.values))

}
