
bs.M.CV <- function(t.grid, z.grid, O.tz.arr, E.tz.arr, bs.grid, hi.tz.par = NULL) {

  M <- length(t.grid)
  N <- length(z.grid)
  L <- dim(O.tz.arr)[3]

  if (!is.null(hi.tz.par) && !is.matrix(hi.tz.par)) {
    hi.tz.par <- matrix(hi.tz.par, nrow = M, ncol = N, byrow = FALSE)
  }

  # Univariate kernel: sextic kernel
  K<-function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1))}
  Tt<-matrix(rep(t.grid, times=M),nrow = M, ncol = M,byrow=FALSE)
  # Tt is an MxM matrix with the values of the grid
  Tx<-Tt-t(Tt)
  Zt <- matrix(rep(z.grid, times=N),nrow = N, ncol = N,byrow=FALSE)
  # Zt is an NxN matrix with the values of the grid
  Zx<-Zt-t(Zt)

  cv.score <- function(b.ind) {

    ## Compute the usual multiplicative bias-corrected hazard estimator with
    ## possible bandiwthd
    bs <- as.numeric(bs.grid[b.ind, ])

    K.bt <- K(Tx/bs[1])/bs[1]
    K.bz <- K(Zx/bs[2])/bs[2]
    K.diag   <- matrix(K.bt[1, 1] * K.bz[1, 1], nrow = M, ncol = N)

    aux.list   <- vector("list", L)
    pilot.list <- vector("list", L)
    alpha.curr <- matrix(1, nrow = M, ncol = N)

    for (k in seq_len(L)) {
      ## Moment calculations involved in the estimator
      pilot.list[[k]] <- alpha.curr
      W <- alpha.curr * E.tz.arr[, , k]

      c0  <- K.bt %*% W %*% t(K.bz)
      c10 <- (K.bt*Tx) %*% W %*% t(K.bz)
      c11 <- K.bt %*% W %*% t((K.bz * Zx))
      d00 <- (K.bt*(Tx^2)) %*% W %*% t(K.bz)
      d01 <- (K.bt*Tx) %*% W %*% t((K.bz * Zx))
      d11 <- K.bt %*% W %*% t((K.bz*(Zx^2)))
      det.D <- d00 * d11 - d01^2
      den   <- det.D * c0 - (c10 * (d11 * c10 - d01 * c11) + c11 * (d00 * c11 - d01 * c10))

      num1 <- (K.bt %*% O.tz.arr[,,k] %*% t(K.bz)) * det.D
      num2 <- ((K.bt*Tx) %*% O.tz.arr[,,k] %*% t(K.bz)) * (d11 * c10 - d01 * c11)
      num3 <- (K.bt %*% O.tz.arr[,,k] %*% t((K.bz * Zx))) * (d00 * c11 - d01 * c10)
      num  <- num1 - num2 - num3
      estim <- num/den
      estim[den == 0 | estim < 0] <- NA

      aux.list[[k]] <- list(num = num, den = den, det.D = det.D)

      if (k < L) {
        alpha.curr <- alpha.curr * estim
        alpha.curr[is.na(alpha.curr)] <- 0
      } else {
        estim.tz.bs <- if (!is.null(hi.tz.par)) alpha.curr * hi.tz.par * estim else alpha.curr * estim
      }
    }

    ## Compute the leave-one-out (LOO) multiplicative bias-corrected hazard estimator
    ind.occ <- (O.tz.arr[,,L] > 0) * 1
    alpha.loo <- matrix(1, nrow = M, ncol = N)

    for (k in seq_len(L)) {
      aux.k <- aux.list[[k]]

      if (k == 1) {
        den.loo <- aux.k$den
      } else {
        delta.W <- (alpha.loo - pilot.list[[k]]) * E.tz.arr[, , k]
        den.loo <- aux.k$den + aux.k$det.D * (K.diag * delta.W)
      }

      num.loo <- aux.k$num - ind.occ * K.diag * aux.k$det.D
      estim.k <- num.loo / den.loo

      valid <- !is.na(estim.k) & (den.loo != 0) & (estim.k >= 0) & !is.na(alpha.loo)
      alpha.loo <- ifelse(valid, alpha.loo * estim.k, NA_real_)
    }

    estim.loo.tz.bs <- if (!is.null(hi.tz.par)) alpha.loo * hi.tz.par else alpha.loo

    ## CV Score
    vec.E.tz <- as.double(E.tz.arr[,,L])
    vec.O.tz <- as.double(O.tz.arr[,,L])
    dif1 <- sum((estim.tz.bs^2) * vec.E.tz, na.rm = TRUE)
    dif2 <- sum(estim.loo.tz.bs * vec.O.tz, na.rm = TRUE)
    cv.b <- dif1 - 2 * dif2
    if (is.infinite(cv.b) || cv.b == 0) cv.b <- NA
    return(cv.b)
  }

  nb <- nrow(bs.grid)
  cv.values <- sapply(seq_len(nb), cv.score)
  ind.cv <- which.min(cv.values)
  if ((ind.cv == 1) || (ind.cv == nb)) {
    warning("The CV score does not have a minimum in the grid of bandwidths")
  }

  return(list(bcv = bs.grid[ind.cv, ], cv.values = cv.values))
}
