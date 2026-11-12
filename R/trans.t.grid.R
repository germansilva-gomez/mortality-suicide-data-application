.loglik <- function(par, t, O, E, n.int = 300){

  a0 <- par[1]
  a1 <- par[2]
  a2 <- par[3]
  sigma <- par[4]

  alpha <- numeric(length(t))

  for(i in seq_along(t)){

    s.grid <- seq(0, t[i], length.out = n.int)
    ds <- diff(s.grid)[1]

    f.vals <- exp(a0 + a1*s.grid + a2*s.grid^2)

    integral.approx <- sum(diff(s.grid) * (head(f.vals,-1) + tail(f.vals,-1))/2)

    alpha[i] <- exp(a0 + a1*t[i] + a2*t[i]^2) / (1 + sigma*integral.approx)
  }

  ll <- sum(O * log(alpha) - E * alpha)

  return(-ll)
}

trans.t.grid <- function(t.grid, O.t, E.t, E.tz, par.ini) {
  # O.t es el vector marginal del tiempo de ocurrencias
  # E.t es el vector marginal del tiempo de expuestos a las ocurrencias

  nrow <- length(t.grid)

  t.grid.modif <- rep(0, nrow)
  hazard.par <- rep(0, nrow)

  res <- solnp(par.ini, .loglik,
               t = t.grid, O = O.t, E = E.t, LB = c(-Inf, 0, 0, 0), UB = c(Inf, Inf, Inf, Inf))

  a0 <- res$par[1]
  a1 <- res$par[2]
  a2 <- res$par[3]
  sigma <- res$par[4]

  alpha.par <- function(t){

    n.int <- 300

    s.grid <- seq(0, t, length.out = n.int)
    ds <- diff(s.grid)[1]

    f.vals <- exp(a0 + a1*s.grid + a2*s.grid^2)

    integral.approx <- sum(diff(s.grid) * (head(f.vals,-1) + tail(f.vals,-1))/2)

    return(exp(a0 + a1*t + a2*t^2) / (1 + sigma*integral.approx))
  }

  alpha.par.eval <- sapply(t.grid, alpha.par)

  hazard.par <- alpha.par.eval

  t.grid.modif <- cumsum(alpha.par.eval)*diff(t.grid)[1]

  t.grid.modif.diff <- diff(c(0,t.grid.modif))

  E.tz.modif <- E.tz*t.grid.modif.diff/diff(t.grid)[1]

  return(list(t.grid.modif=t.grid.modif, E.tz.modif=E.tz.modif, hi.tz.par = hazard.par, res = res))
}
