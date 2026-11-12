
b.M.CV.4Levels <- function(t.grid, z.grid, O.tz, O.tz.1, O.tz.2, O.tz.3,
                                       E.tz, E.tz.1, E.tz.2, E.tz.3, bs.grid,
                                       trans = FALSE, hi.tz.par = NULL)
{
  M <- length(t.grid)
  N <- length(z.grid)

  if (!is.matrix(hi.tz.par) && trans == TRUE) {
    hi.tz.par <- matrix(hi.tz.par, nrow = M, ncol = N, byrow = FALSE)
  }

  hazard.b <- function(t.grid, z.grid, O.tz, E.tz, alpha.x, alpha.grid, bs, ...) {

    # Univariate kernel: sextic kernel
    K <- function(u) { return(((3003/2048)*(1-(u)^2)^6)*(abs(u)<1)) }

    ## 1. Compute kernel evaluations and moments of the kernel
    ## 1.1. Kernel evaluations at dimension t
    Tt <- matrix(rep(t.grid, times=M), nrow = M, ncol = M, byrow=FALSE)
    Tx <- Tt - t(Tt)
    ut <- Tx / bs[1,1]
    K.bt <- apply(ut, 1:2, K) / bs[1,1]

    ## 1.2. Dimension z
    Zt <- matrix(rep(z.grid, times=N), nrow = N, ncol = N, byrow=FALSE)
    Zx <- Zt - t(Zt)
    ut <- Zx / bs[1,2]
    K.bz <- apply(ut, 1:2, K) / bs[1,2]
    rm(ut, Tt, Zt)

    c0  <- (K.bt) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz)
    c10 <- (K.bt * Tx) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz)
    c11 <- (K.bt) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz*Zx)
    d00 <- (K.bt * (Tx^2)) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz)
    d01 <- (K.bt * Tx) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz*Zx)
    d11 <- (K.bt) %*% ((alpha.grid) * (E.tz)) %*% t(K.bz*(Zx^2))

    det.D <- d00*d11 - d01^2
    den <- det.D*c0 - ( c10*(d11*c10 - d01*c11) + c11*(d00*c11 - d01*c10) )

    num1 <- ((K.bt) %*% (O.tz) %*% t(K.bz)) * (det.D)
    num2 <- ((K.bt * Tx) %*% (O.tz) %*% t(K.bz)) * (d11*c10 - d01*c11)
    num3 <- ((K.bt) %*% (O.tz) %*% t(K.bz*Zx)) * (d00*c11 - d01*c10)

    num <- num1 - num2 - num3
    estim <- num / den
    estim[den == 0] <- NA
    estim[estim < 0] <- NA

    g.tz <- matrix(as.double(estim), M, N, byrow=FALSE)
    g.tz[which(g.tz < 0)] <- NA

    hi.tz <- alpha.x * g.tz

    if (!is.null(hi.tz.par)) {
      num.trans <- alpha.x * hi.tz.par * num
      estim.trans <- num.trans / den
      estim.trans[den == 0] <- NA
      estim.trans[estim.trans < 0] <- NA
      return(list(K.bt = K.bt, K.bz = K.bz, Tx = Tx, Zx = Zx, c0 = c0, c10 = c10,
                  c11 = c11, d00 = d00, d01 = d01, d11 = d11, det.D = det.D,
                  den = den, hi.tz = estim.trans, hi.tz.trans = hi.tz))
    } else {
      return(list(K.bt = K.bt, K.bz = K.bz, Tx = Tx, Zx = Zx, c0 = c0, c10 = c10,
                  c11 = c11, d00 = d00, d01 = d01, d11 = d11, det.D = det.D,
                  den = den, hi.tz = hi.tz))
    }
  }

  ## CV.score
  cv.score <- function(b.ind) {

    if (is.null(hi.tz.par)) {
    ## The usual estimator (not loo) for the nonparametric hazard estimator
    # Level 0
    aux.0 <- hazard.b(t.grid, z.grid, O.tz, E.tz, matrix(1, M, N), matrix(1, M, N), bs.grid[b.ind, ])
    alpha.pilot.0 <- aux.0$hi.tz
    alpha.pilot.0[is.na(alpha.pilot.0)] <- 0

    # Level 1
    aux.1 <- hazard.b(t.grid, z.grid, O.tz.1, E.tz.1, alpha.pilot.0, alpha.pilot.0, bs.grid[b.ind, ])
    alpha.pilot.1 <- aux.1$hi.tz
    alpha.pilot.1[is.na(alpha.pilot.1)] <- 0

    # Level 2
    aux.2 <- hazard.b(t.grid, z.grid, O.tz.2, E.tz.2, alpha.pilot.1, alpha.pilot.1, bs.grid[b.ind, ])
    alpha.pilot.2 <- aux.2$hi.tz
    alpha.pilot.2[is.na(alpha.pilot.2)] <- 0

    # Level 3
    aux.3 <- hazard.b(t.grid, z.grid, O.tz.3, E.tz.3, alpha.pilot.2, alpha.pilot.2, bs.grid[b.ind, ])
    estim.tz.bs <- aux.3$hi.tz

    } else {

    ## The usual estimator (not loo) for the semiparametric hazard estimator
    # Level 0
    aux.0 <- hazard.b(t.grid, z.grid, O.tz, E.tz, matrix(1, M, N), matrix(1, M, N), bs.grid[b.ind, ], hi.tz.par)
    alpha.pilot.0 <- aux.0$hi.tz.trans
    alpha.pilot.0[is.na(alpha.pilot.0)] <- 0

    # Level 1
    aux.1 <- hazard.b(t.grid, z.grid, O.tz.1, E.tz.1, alpha.pilot.0, alpha.pilot.0, bs.grid[b.ind, ], hi.tz.par)
    alpha.pilot.1 <- aux.1$hi.tz.trans
    alpha.pilot.1[is.na(alpha.pilot.1)] <- 0

    # Level 2
    aux.2 <- hazard.b(t.grid, z.grid, O.tz.2, E.tz.2, alpha.pilot.1, alpha.pilot.1, bs.grid[b.ind, ], hi.tz.par)
    alpha.pilot.2 <- aux.2$hi.tz.trans
    alpha.pilot.2[is.na(alpha.pilot.2)] <- 0

    # Level 3
    aux.3 <- hazard.b(t.grid, z.grid, O.tz.3, E.tz.3, alpha.pilot.2, alpha.pilot.2, bs.grid[b.ind, ], hi.tz.par)
    estim.tz.bs <- aux.3$hi.tz

    }

    calc.moments.ij <- function(i, j, K.bt, K.bz, Tx, Zx, alpha.grid, E.tz) {
      W <- alpha.grid * E.tz
      Kb.t_i <- K.bt[i, ]
      Kb.z_j <- K.bz[j, ]
      Tx_i   <- Tx[i, ]
      Zx_j   <- Zx[j, ]

      c0  <- as.numeric(Kb.t_i %*% W %*% Kb.z_j)
      c10 <- as.numeric((Kb.t_i * Tx_i) %*% W %*% Kb.z_j)
      c11 <- as.numeric(Kb.t_i %*% W %*% (Kb.z_j * Zx_j))
      d00 <- as.numeric((Kb.t_i * (Tx_i^2)) %*% W %*% Kb.z_j)
      d01 <- as.numeric((Kb.t_i * Tx_i) %*% W %*% (Kb.z_j * Zx_j))
      d11 <- as.numeric(Kb.t_i %*% W %*% (Kb.z_j * (Zx_j^2)))

      det.D <- d00 * d11 - d01^2
      den   <- det.D * c0 - (c10 * (d11 * c10 - d01 * c11) + c11 * (d00 * c11 - d01 * c10))

      return(list(c0 = c0, c10 = c10, c11 = c11,
                  d00 = d00, d01 = d01, d11 = d11,
                  det.D = det.D, den = den))
    }

    ## The leave-one-out (loo) hazard estimator
    ## Create a function similar to hazard.b but leave-one-out
    hazard.loo.4levels <- function(aux.0, aux.1, aux.2, aux.3,
                                   O.tz, O.tz.1, O.tz.2, O.tz.3,
                                   E.tz.1, E.tz.2, E.tz.3,
                                   alpha.pilot.0, alpha.pilot.1, alpha.pilot.2,
                                   M, N, ...)
    {
      estim.ij <- double(M * N)
      l <- 0

      alpha.grid.1 <- alpha.pilot.0
      alpha.grid.2 <- alpha.pilot.1
      alpha.grid.3 <- alpha.pilot.2

      for (i in 1:M)
      {
        for (j in 1:N)
        {
          l <- l + 1
          ij <- 0
          if (O.tz.3[i,j] > 0) {
            O.tz.3[i,j] <- O.tz.3[i,j] - 1
            O.tz.2[i,j] <- O.tz.2[i,j] - 1
            O.tz.1[i,j] <- O.tz.1[i,j] - 1
            O.tz[i,j]   <- O.tz[i,j]   - 1
            ij <- 1
          }

          # alpha_0^{-i}(x_i)
          num1.0 <- ((aux.0$K.bt[i,]) %*% (O.tz) %*% (aux.0$K.bz[j,])) * (aux.0$det.D[i,j])
          num2.0 <- ((aux.0$K.bt[i,] * aux.0$Tx[i,]) %*% (O.tz) %*% aux.0$K.bz[j,]) * (aux.0$d11[i,j]*aux.0$c10[i,j] - aux.0$d01[i,j]*aux.0$c11[i,j])
          num3.0 <- ((aux.0$K.bt[i,]) %*% (O.tz) %*% (aux.0$K.bz[j,]*aux.0$Zx[j,])) * (aux.0$d00[i,j]*aux.0$c11[i,j] - aux.0$d01[i,j]*aux.0$c10[i,j])
          num.0  <- num1.0 - num2.0 - num3.0
          estim.0 <- num.0 / aux.0$den[i,j]

          if (!is.na(estim.0) && aux.0$den[i,j] != 0 && estim.0 >= 0) {
            alpha.0.loo <- as.numeric(estim.0)
          } else {
            alpha.0.loo <- NA
          }

          # alpha_1^{-i}(x_i)
          if (!is.na(alpha.0.loo)) {
            alpha.grid.1[i,j] <- alpha.0.loo

            mom.1 <- calc.moments.ij(i, j, aux.1$K.bt, aux.1$K.bz, aux.1$Tx, aux.1$Zx, alpha.grid.1, E.tz.1)

            num1.1 <- ((aux.1$K.bt[i,]) %*% (O.tz.1) %*% (aux.1$K.bz[j,])) * (mom.1$det.D)
            num2.1 <- ((aux.1$K.bt[i,] * aux.1$Tx[i,]) %*% (O.tz.1) %*% aux.1$K.bz[j,]) * (mom.1$d11 * mom.1$c10 - mom.1$d01 * mom.1$c11)
            num3.1 <- ((aux.1$K.bt[i,]) %*% (O.tz.1) %*% (aux.1$K.bz[j,]*aux.1$Zx[j,])) * (mom.1$d00 * mom.1$c11 - mom.1$d01 * mom.1$c10)
            num.1  <- num1.1 - num2.1 - num3.1
            estim.1 <- num.1 / mom.1$den

            if (!is.na(estim.1) && mom.1$den != 0 && estim.1 >= 0) {
              alpha.1.loo <- alpha.0.loo * as.numeric(estim.1)
            } else {
              alpha.1.loo <- NA
            }

            alpha.grid.1[i,j] <- alpha.pilot.0[i,j]
          } else {
            alpha.1.loo <- NA
          }

          # alpha_2^{-i}(x_i)
          if (!is.na(alpha.1.loo)) {
            alpha.grid.2[i,j] <- alpha.1.loo

            mom.2 <- calc.moments.ij(i, j, aux.2$K.bt, aux.2$K.bz, aux.2$Tx, aux.2$Zx, alpha.grid.2, E.tz.2)

            num1.2 <- ((aux.2$K.bt[i,]) %*% (O.tz.2) %*% (aux.2$K.bz[j,])) * (mom.2$det.D)
            num2.2 <- ((aux.2$K.bt[i,] * aux.2$Tx[i,]) %*% (O.tz.2) %*% aux.2$K.bz[j,]) * (mom.2$d11 * mom.2$c10 - mom.2$d01 * mom.2$c11)
            num3.2 <- ((aux.2$K.bt[i,]) %*% (O.tz.2) %*% (aux.2$K.bz[j,]*aux.2$Zx[j,])) * (mom.2$d00 * mom.2$c11 - mom.2$d01 * mom.2$c10)
            num.2  <- num1.2 - num2.2 - num3.2
            estim.2 <- num.2 / mom.2$den

            if (!is.na(estim.2) && mom.2$den != 0 && estim.2 >= 0) {
              alpha.2.loo <- alpha.1.loo * as.numeric(estim.2)
            } else {
              alpha.2.loo <- NA
            }

            alpha.grid.2[i,j] <- alpha.pilot.1[i,j]
          } else {
            alpha.2.loo <- NA
          }

          # alpha_3^{-i}(x_i)
          if (!is.na(alpha.2.loo)) {
            alpha.grid.3[i,j] <- alpha.2.loo

            mom.3 <- calc.moments.ij(i, j, aux.3$K.bt, aux.3$K.bz, aux.3$Tx, aux.3$Zx, alpha.grid.3, E.tz.3)

            num1.3 <- ((aux.3$K.bt[i,]) %*% (O.tz.3) %*% (aux.3$K.bz[j,])) * (mom.3$det.D)
            num2.3 <- ((aux.3$K.bt[i,] * aux.3$Tx[i,]) %*% (O.tz.3) %*% aux.3$K.bz[j,]) * (mom.3$d11 * mom.3$c10 - mom.3$d01 * mom.3$c11)
            num3.3 <- ((aux.3$K.bt[i,]) %*% (O.tz.3) %*% (aux.3$K.bz[j,]*aux.3$Zx[j,])) * (mom.3$d00 * mom.3$c11 - mom.3$d01 * mom.3$c10)
            num.3  <- num1.3 - num2.3 - num3.3
            estim.3 <- num.3 / mom.3$den

            if (!is.na(estim.3) && mom.3$den != 0 && estim.3 >= 0) {
              if (is.null(hi.tz.par)) {estim.ij[l] <- alpha.2.loo * as.numeric(estim.3) }
              else {estim.ij[l] <- alpha.2.loo * as.numeric(estim.3) * hi.tz.par[i,j]}
            } else {
              estim.ij[l] <- NA
            }

            alpha.grid.3[i,j] <- alpha.pilot.2[i,j]
          } else {
            estim.ij[l] <- NA
          }

          if (ij > 0) {
            O.tz.3[i,j] <- O.tz.3[i,j] + 1
            O.tz.2[i,j] <- O.tz.2[i,j] + 1
            O.tz.1[i,j] <- O.tz.1[i,j] + 1
            O.tz[i,j]   <- O.tz[i,j]   + 1
          }
        }
      }
      return(estim.ij)
    }

    if (!is.null(hi.tz.par)) {

      loo_vec <- hazard.loo.4levels(aux.0, aux.1, aux.2, aux.3,
                                    O.tz, O.tz.1, O.tz.2, O.tz.3,
                                    E.tz.1, E.tz.2, E.tz.3,
                                    alpha.pilot.0, alpha.pilot.1,
                                    alpha.pilot.2,
                                    M, N, hi.tz.par)

    } else {

      loo_vec <- hazard.loo.4levels(aux.0, aux.1, aux.2, aux.3,
                                    O.tz, O.tz.1, O.tz.2, O.tz.3,
                                    E.tz.1, E.tz.2, E.tz.3,
                                    alpha.pilot.0, alpha.pilot.1,
                                    alpha.pilot.2,
                                    M, N)
    }

    estim.loo.tz.bs <- matrix(loo_vec, nrow = M, ncol = N, byrow = TRUE)

    # Score CV
    vec.E.tz <- as.double(E.tz.3)
    vec.O.tz <- as.double(O.tz.3)
    dif1 <- sum((estim.tz.bs^2) * vec.E.tz, na.rm = TRUE)
    dif2 <- sum(estim.loo.tz.bs * vec.O.tz, na.rm = TRUE)
    cv.b <- dif1 - 2 * dif2
    if (is.infinite(cv.b) || cv.b == 0) cv.b <- NA

    return(cv.b)
  }

  nb <- nrow(bs.grid)
  cv.values <- sapply(1:nb, cv.score)
  ind.cv <- which.min(cv.values)
  if ((ind.cv == 1) || (ind.cv == nb))
    warning("The CV score does not have a minimum in the grid of bandwidths")
  bcv <- bs.grid[ind.cv, ]

  return(list(bcv = bcv, cv.values = cv.values))
}
