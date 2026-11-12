hazard2D.M<-function(t.grid,z.grid,O.tz,E.tz,alpha.grid,alpha.x,bs,trans = FALSE,hi.tz.par = NULL)

{
  # Univariate kernel: sextic kernel
  K<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}

  ## 1. Compute kernel evaluations and moments of the kernel
  ## 1.1. Kernel evaluations at dimension t
  M<-length(t.grid)
  Tt<-matrix(rep(t.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
  # Tt is an MxM matrix with the values of the grid
  Tx<-Tt-t(Tt)
  nb<-nrow(bs)
  K.bt<-array(0,dim=c(M,M,nb))
  # Each element of K.bt is the Kernel evaluated at ut below
  ut<-array(0,dim=c(M,M,nb))
  for(b in 1:nb)
  {
    ut[,,b]=Tx/bs[b,1]
    K.bt[,,b]<-apply(ut[,,b],1:2,K)
    # This works and returns the value of the Kernel function at each individual value of the array ut
    K.bt[,,b]<-K.bt[,,b]/(bs[b,1])
  }

  ## 1.2. Dimension z (evaluation points are the same but bandwidth not necessarily)
  N <- length(z.grid)
  Zt <- matrix(rep(z.grid, times=N),nrow = N, ncol = N,byrow=FALSE)
  # Zt is an NxN matrix with the values of the grid
  Zx<-Zt-t(Zt)
  K.bz<-array(0,dim=c(N,N,nb))
  # Each element of K.bz is the Kernel evaluated at ut below
  uz<-array(0,dim=c(N,N,nb))
  for(b in 1:nb)
  {
    uz[,,b]=Zx/bs[b,2]
    K.bz[,,b]<-apply(uz[,,b],1:2,K)
    # This works and returns the value of the Kernel function at each individual value of the array ut
    K.bz[,,b]<-K.bz[,,b]/(bs[b,2])
  }
  rm(ut,uz,Tt,Zt)

  ## 1.3. Moment calculations involved in the estimator
  c0<-c10<-c11<-d00<-d01<-d11<-det.D<-den<-array(0,dim=c(M,N,nb))
  for (b in 1:nb)
  {

    c0[,,b]<-(K.bt[,,b]) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz[,,b])
    c10[,,b]<-(K.bt[,,b] * Tx) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz[,,b])
    c11[,,b]<-(K.bt[,,b]) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz[,,b]*Zx)
    d00[,,b]<-(K.bt[,,b] * (Tx^2)) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz[,,b])
    d01[,,b]<-(K.bt[,,b] * Tx) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz[,,b]*Zx)
    d11[,,b]<-(K.bt[,,b]) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz[,,b]*(Zx^2))

  }

  ## 2. Create a function to compute the hazard from results above
  hazard.b <- function(K.bt, K.bz, Tx, Zx, b, O.tz, c10, c11, d00, d01, d11, det.D, den, ...) {

    det.D[,,b] <- d00[,,b]*d11[,,b] - d01[,,b]^2
    den[,,b]   <- det.D[,,b]*c0[,,b] -
      ( c10[,,b]*(d11[,,b]*c10[,,b] - d01[,,b]*c11[,,b]) +
          c11[,,b]*(d00[,,b]*c11[,,b] - d01[,,b]*c10[,,b]) )

    num1 <- ((K.bt[,,b]) %*% (O.tz) %*% t(K.bz[,,b])) * (det.D[,,b])
    num2 <- ((K.bt[,,b] * Tx) %*% (O.tz) %*% t(K.bz[,,b])) * (d11[,,b]*c10[,,b] - d01[,,b]*c11[,,b])
    num3 <- ((K.bt[,,b]) %*% (O.tz) %*% t(K.bz[,,b]*Zx)) * (d00[,,b]*c11[,,b] - d01[,,b]*c10[,,b])

    num <- num1 - num2 - num3
    estim <- num / den[,,b]
    estim[den[,,b] == 0] <- NA
    estim[estim < 0]     <- NA

    if (is.null(hi.tz.par) && isTRUE(trans)) {stop("argument 'hi.tz.par' is required when 'trans = TRUE'")}
    if (!is.null(hi.tz.par) && isFALSE(trans)) {
      stop("Conflict: 'hi.tz.par' provided but 'trans' is FALSE")
    }

    if (!is.null(hi.tz.par)) {
      num.trans <- alpha.x * hi.tz.par * num
      estim.trans <- num.trans / den[,,b]
      estim.trans[den[,,b] == 0] <- NA
      estim.trans[estim.trans < 0] <- NA
      return(list(estim = estim.trans, estim.trans = estim))
    } else {
      return(estim)
    }
  }

  ## 3. Calculate the nonparametric or the semiparametric version of the estimator
  if (is.null(hi.tz.par) && isFALSE(trans)) {

    estim.tz.b <- hazard.b(K.bt, K.bz, Tx, Zx, b=1, O.tz, c10, c11, d00, d01, d11, det.D, den)
    g.tz <- matrix(as.double(estim.tz.b), M, N, byrow=FALSE)
    g.tz[which(g.tz < 0)] <- NA
    hi.tz <- alpha.x * g.tz

    return(list(hi.tz = hi.tz))

  } else {

    estim.tz.b <- hazard.b(K.bt, K.bz, Tx, Zx, b=1, O.tz, c10, c11, d00, d01, d11, det.D, den,
                           alpha.x, hi.tz.par)

    g.tz <- matrix(as.double(estim.tz.b$estim.trans), M, N, byrow=FALSE)
    g.tz[which(g.tz < 0)] <- NA
    hi.tz.trans <- alpha.x * g.tz

    hi.tz <- matrix(as.double(estim.tz.b$estim), M, N, byrow=FALSE)

    return(list(hi.tz = hi.tz, hi.tz.trans = hi.tz.trans))
  }
}
