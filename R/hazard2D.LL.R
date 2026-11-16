hazard2D.LL <- function(t.grid, z.grid, O.tz, E.tz, bs.grid, cv = TRUE) {

  M  <- length(t.grid)
  N  <- length(z.grid)
  nb <- nrow(bs.grid)

  # Univariate kernel: sextic kernel
  K<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}

  ## 1. Compute kernel evaluations and moments of the kernel
  ## 1.1. Kernel evaluations at dimension t
  Tt<-matrix(rep(t.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
  # Tt is an MxM matrix with the values of the grid
  Tx<-Tt-t(Tt)
  Zt <- matrix(rep(z.grid, times=N),nrow = N, ncol = N,byrow=FALSE)
  # Zt is an NxN matrix with the values of the grid
  Zx<-Zt-t(Zt)

  ## 1.2. Moment calculations involved in the estimator
  hazard.b <- function(b) {
    bs <- as.numeric(bs.grid[b, ])
    K.bt <- K(Tx/bs[1])/bs[1]
    K.bz <- K(Zx/bs[2])/bs[2]

    c0  <- K.bt %*% E.tz %*% t(K.bz)
    c10 <- (K.bt*Tx) %*% E.tz %*% t(K.bz)
    c11 <- K.bt %*% E.tz %*% t(K.bz*Zx)
    d00 <- (K.bt*Tx^2) %*% E.tz %*% t(K.bz)
    d01 <- (K.bt*Tx) %*% E.tz %*% t(K.bz*Zx)
    d11 <- K.bt %*% E.tz %*% t(K.bz * Zx^2)

    det.D <- d00 *d11  - d01 ^2
    den   <- det.D*c0 - (c10*(d11*c10 - d01*c11) + c11*(d00*c11 - d01*c10) )

    num1 <- ((K.bt) %*% (O.tz) %*% t(K.bz)) * (det.D)
    num2 <- ((K.bt * Tx) %*% (O.tz) %*% t(K.bz)) * (d11*c10 - d01*c11)
    num3 <- ((K.bt) %*% (O.tz) %*% t(K.bz*Zx)) * (d00*c11 - d01*c10)
    num <- num1 - num2 - num3
    estim <- num/den
    estim[den == 0] <- NA
    estim[estim < 0] <- NA

    return(list(estim = estim, num = num, den = den, det.D  = det.D,
                K.diag = K.bt[1, 1] * K.bz[1, 1]))
  }

  cv.values <- NULL

  if (cv == TRUE && nb > 1) {

    ## Compute the leave-one-out (LOO) multiplicative bias-corrected hazard estimator
    ind.occ  <- (O.tz > 0) * 1
    vec.E.tz <- as.double(E.tz)
    vec.O.tz <- as.double(O.tz)

    estim.list <- vector("list", nb)
    cv.values  <- numeric(nb)

    for (b in seq_len(nb)) {
      res.b <- hazard.b(b)
      estim.list[[b]] <- res.b$estim

      num.loo   <- res.b$num - ind.occ * res.b$K.diag * res.b$det.D
      estim.loo <- num.loo/res.b$den
      estim.loo[res.b$den == 0] <- NA

      dif1.b <- sum((as.double(res.b$estim))^2 * vec.E.tz, na.rm = TRUE)
      dif2.b <- sum(as.double(estim.loo) * vec.O.tz, na.rm = TRUE)
      cv.b   <- dif1.b- 2*dif2.b

      if (is.infinite(cv.b) || cv.b == 0) cv.b <- NA
      cv.values[b] <- cv.b
    }

    bcv <- which.min(cv.values)
    hi.tz <- estim.list[[bcv]]

  } else {

    bcv   <- 1
    hi.tz <- hazard.b(b = 1)$estim

  }

  return(list(hi.tz = hi.tz, bcv = bs.grid[bcv, ], cv.values = cv.values))
}
