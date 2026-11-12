bs.M.CV.2Levels <- function(t.grid, z.grid, O.tz, O.tz.1, E.tz, E.tz.1, bs.grid, trans = FALSE, hi.tz.par = NULL)
{

  if (is.null(hi.tz.par) && isTRUE(trans)) {stop("argument 'hi.tz.par' is required when 'trans = TRUE'")}
  if (!is.null(hi.tz.par) && isFALSE(trans)) {
    stop("Conflict: 'hi.tz.par' provided but 'trans' is FALSE")
  }

  M<-length(t.grid)
  N <- length(z.grid)

  hazard.b <- function(t.grid, z.grid, O.tz, E.tz, alpha.x, alpha.grid, bs, ...) {

    # Univariate kernel: sextic kernel
    K<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}

    ## 1. Compute kernel evaluations and moments of the kernel
    ## 1.1. Kernel evaluations at dimension t
    Tt<-matrix(rep(t.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
    # Tt is an MxM matrix with the values of the grid
    Tx<-Tt-t(Tt)
    K.bt<-matrix(0,M,M)
    # each element of K.bt is the Kernel evaluated at ut below
    ut<-matrix(0,M,M)
    ut=Tx/(bs[1,1])
    K.bt<-apply(ut,1:2,K)
    K.bt<-K.bt/(bs[1,1])

    ## 1.2. Dimension z
    Zt <- matrix(rep(z.grid, times=N),nrow = N, ncol = N,byrow=FALSE)
    # Zt is an NxN matrix with the values of the grid
    Zx<-Zt-t(Zt)
    K.bz<-matrix(0,N,N)
    ut<-matrix(0,N,N)
    ut=Zx/(bs[1,2])
    K.bz<-apply(ut,1:2,K)
    K.bz<-K.bz/(bs[1,2])
    rm(ut,Tt,Zt)

    c0<-c10<-c11<-d00<-d01<-d11<-det.D<-den<-matrix(0,M,N)

    c0<-(K.bt) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz)
    c10<-(K.bt * Tx) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz)
    c11<-(K.bt) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz*Zx)
    d00<-(K.bt * (Tx^2)) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz)
    d01<-(K.bt * Tx) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz*Zx)
    d11<-(K.bt) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz*(Zx^2))

    det.D<- d00*d11-d01^2
    den<- det.D*c0 -
      ( c10*(d11*c10-d01*c11)
        + c11*(d00*c11 - d01*c10) )

    num1<-((K.bt) %*% (O.tz)  %*%  t(K.bz)) *(det.D)
    num2<-( (K.bt * Tx) %*% (O.tz)  %*%  t(K.bz))*(d11*c10 - d01*c11)
    num3<-( (K.bt) %*% (O.tz)  %*%  t(K.bz*Zx))*(d00*c11 - d01*c10)

    num<-num1-num2-num3
    estim<-num/den
    estim[den==0]<-NA
    estim[estim<0]<-NA

    g.tz<-matrix(as.double(estim),M,N,byrow=FALSE)
    g.tz[which(g.tz<0)]<-NA

    hi.tz <- alpha.x * g.tz

    if (!is.null(hi.tz.par)) {
      num.trans <- alpha.x * hi.tz.par * num
      estim.trans <- num.trans / den
      estim.trans[den == 0] <- NA
      estim.trans[estim.trans < 0] <- NA
      return(list(K.bt = K.bt, K.bz = K.bz, Tx = Tx, Zx = Zx, c0 = c0, c10 = c10,
                  c11 = c11, d00 = d00, d01 = d01, d11 = d11, det.D = det.D,
                  den = den, hi.tz = estim.trans))
    } else {
      return(list(K.bt = K.bt, K.bz = K.bz, Tx = Tx, Zx = Zx, c0 = c0, c10 = c10,
                  c11 = c11, d00 = d00, d01 = d01, d11 = d11, det.D = det.D,
                  den = den, hi.tz = hi.tz))
    }
  }

  ## CV.score
  cv.score <- function(b.ind) {

    ## The leave-one-out (loo) hazard estimator
    ## Create a function similar to hazard.b but leave-one-out
    hazard.loo.b<-function(K.bt,K.bz,Tx,Zx,O.tz,c10,c11,d00,d01,d11,
                           det.D,den,M,N,alpha.x)
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
          num1<-( (K.bt[i,]) %*% (O.tz)  %*%  (K.bz[j,])) *(det.D[i,j])
          num2<-( (K.bt[i,] * Tx[i,]) %*% (O.tz)  %*%  K.bz[j,])*(d11[i,j]*c10[i,j] - d01[i,j]*c11[i,j])
          num3<-( (K.bt[i,]) %*% (O.tz)  %*%  (K.bz[j,]*Zx[j,]) )*(d00[i,j]*c11[i,j] - d01[i,j]*c10[i,j])
          num<-num1-num2-num3
          estim<-num/den[i,j]; estim[den[i,j]==0]<-NA
          if (!is.na(estim)) {
            estim.ij[l] <- alpha.x[i,j] * estim
          } else {
            estim.ij[l] <- NA
          }
          if (ij>0) {O.tz[i,j]<-Oij.prev; ij=0}
        }
      }
      return(estim.ij)
    }

    ## Evaluate the function above
    alpha.pilot <- hazard2D.LL(t.grid=t.grid,z.grid=z.grid,O.tz=O.tz,E.tz=E.tz,
                               bs.grid = bs.grid[b.ind, ],cv=F)$hi.tz

    alpha.pilot[is.na(alpha.pilot)] <- 0

    aux <- hazard.b(t.grid, z.grid, O.tz.1, E.tz.1, alpha.pilot, alpha.pilot, bs.grid[b.ind, ])

    loo_vec <-hazard.loo.b(K.bt=aux$K.bt,K.bz=aux$K.bz,Tx=aux$Tx,Zx=aux$Zx,
                           O.tz.1,c10=aux$c10,c11=aux$c11,d00=aux$d00,
                           d01=aux$d01,d11=aux$d11,det.D=aux$det.D,den=aux$den,M,N,alpha.pilot)

    estim.loo.tz.bs <- matrix(loo_vec, nrow = M, ncol = N, byrow = TRUE)

    ## The usual estimator (not loo)
    estim.tz.bs <- hazard.b(t.grid, z.grid, O.tz.1, E.tz.1, alpha.pilot, alpha.pilot, bs.grid[b.ind, ])$hi.tz


    if (!is.null(hi.tz.par)) {
      aux <- hazard.b(t.grid, z.grid, O.tz.1, E.tz.1, alpha.pilot, alpha.pilot, bs.grid[b.ind, ], hi.tz.par)

      loo_vec <-hazard.loo.b(K.bt=aux$K.bt,K.bz=aux$K.bz,Tx=aux$Tx,Zx=aux$Zx,
                             O.tz.1,c10=aux$c10,c11=aux$c11,d00=aux$d00,
                             d01=aux$d01,d11=aux$d11,det.D=aux$det.D,den=aux$den,M,N,alpha.pilot)

      estim.loo.tz.bs <- matrix(loo_vec, nrow = M, ncol = N, byrow = TRUE)
      estim.loo.tz.bs.par <- estim.loo.tz.bs * hi.tz.par

      ## The usual estimator (not loo)
      estim.tz.bs <- hazard.b(t.grid, z.grid, O.tz.1, E.tz.1, alpha.pilot, alpha.pilot, bs.grid[b.ind, ], hi.tz.par)$hi.tz
    }

    # CV.b
    vec.E.tz <- E.tz.1
    vec.O.tz <- O.tz.1
    dif1 <- sum((estim.tz.bs^2) * vec.E.tz,na.rm=TRUE)
    dif2 <- sum(estim.loo.tz.bs * (vec.O.tz),na.rm=TRUE)
    cv.b <- dif1 - 2 * dif2
    if (cv.b == Inf | cv.b == 0) cv.b <- NA

    return(cv.b)

  }
  nb <- nrow(bs.grid)
  cv.values <- sapply(1:nb, cv.score)
  ind.cv <- which.min(cv.values)
  if ((ind.cv == 1) | (ind.cv == nb))
    warning("The CV score does not have a minumum in the grid of bandwidths")
  bcv <- bs.grid[ind.cv,]

  return(list(bcv = bcv, cv.values = cv.values))
}
