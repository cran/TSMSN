# Auxiliary functions ---------------------------

# Moments ---------------------------
.mSNI <- function(mu, sigma2, lambda, nu, a, b, type) {

  if (type == "SN") { type1 <- "N" }
  if (type == "ST") { type1 <- "T" }
  if (type == "SSL") { type1 <- "SL" }
  if (type == "SCN") { type1 <- "CN" }

  Lim11 <- (a - mu) / sqrt(sigma2)
  Lim21 <- (b - mu) / sqrt(sigma2)

  a_alpha <- Lim11 * sqrt(1 + lambda^2)
  b_alpha <- Lim21 * sqrt(1 + lambda^2)

  tau <- (.cdfSNI(y = Lim21, mu = 0, sigma2 = 1, lambda = lambda, nu = nu, type = type) -
            .cdfSNI(y = Lim11, mu = 0, sigma2 = 1, lambda = lambda, nu = nu, type = type)) ^ (-1)

  EUX1 <-
    tau * (
      .L(s = 1, lambda = lambda) *
        (.E_Phi(r = -0.5, a = b_alpha, nu = nu, type = type1) -
           .E_Phi(r = -0.5, a = a_alpha, nu = nu, type = type1)) -
        (.E_phiSNI(r = -0.5, a = Lim21, lambda = lambda, nu = nu, type = type) -
           .E_phiSNI(r = -0.5, a = Lim11, lambda = lambda, nu = nu, type = type))
    )

  if (type == "SN" || type == "SCN" || (type == "ST" & nu > 2) || (type == "SSL" & nu > 1)) {
    EUX2 <-
      tau * (
        (.E_PhiSNI(r = -1, a = Lim21, lambda = lambda, nu = nu, type = type) -
           .E_PhiSNI(r = -1, a = Lim11, lambda = lambda, nu = nu, type = type)) -
          .L(s = 2, lambda = lambda) *
          (.E_phi(r = -1, a = b_alpha, nu = nu, type = type1) -
             .E_phi(r = -1, a = a_alpha, nu = nu, type = type1)) -
          (Lim21 * .E_phiSNI(r = -0.5, a = Lim21, lambda = lambda, nu = nu, type = type) -
             Lim11 * .E_phiSNI(r = -0.5, a = Lim11, lambda = lambda, nu = nu, type = type))
      )
  }

  if (type == "SN" || type == "SCN" || (type == "ST" & nu > 3) || (type == "SSL" & nu > 1.5)) {
    EUX3 <-
      tau * (
        2 * .L(s = 1, lambda = lambda) *
          (.E_Phi(r = -1.5, a = b_alpha, nu = nu, type = type1) -
             .E_Phi(r = -1.5, a = a_alpha, nu = nu, type = type1)) -
          2 * (.E_phiSNI(r = -1.5, a = Lim21, lambda = lambda, nu = nu, type = type) -
                 .E_phiSNI(r = -1.5, a = Lim11, lambda = lambda, nu = nu, type = type)) -
          ((Lim21^2) * .E_phiSNI(r = -0.5, a = Lim21, lambda = lambda, nu = nu, type = type) -
             (Lim11^2) * .E_phiSNI(r = -0.5, a = Lim11, lambda = lambda, nu = nu, type = type)) +
          .L(s = 3, lambda = lambda) *
          (.E_Phi(r = -1.5, a = b_alpha, nu = nu, type = type1) -
             .E_Phi(r = -1.5, a = a_alpha, nu = nu, type = type1)) -
          .L(s = 2, lambda = lambda) *
          (Lim21 * .E_phi(r = -1, a = b_alpha, nu = nu, type = type1) -
             Lim11 * .E_phi(r = -1, a = a_alpha, nu = nu, type = type1))
      )
  }

  if (type == "SN" || type == "SCN" || (type == "ST" & nu > 4) || (type == "SSL" & nu > 2)) {
    EUX4 <-
      tau * (
        3 * (.E_PhiSNI(r = -2, a = Lim21, lambda = lambda, nu = nu, type = type) -
               .E_PhiSNI(r = -2, a = Lim11, lambda = lambda, nu = nu, type = type)) -
          3 * .L(s = 2, lambda = lambda) *
          (.E_phi(r = -2, a = b_alpha, nu = nu, type = type1) -
             .E_phi(r = -2, a = a_alpha, nu = nu, type = type1)) -
          3 * (Lim21 * .E_phiSNI(r = -1.5, a = Lim21, lambda = lambda, nu = nu, type = type) -
                 Lim11 * .E_phiSNI(r = -1.5, a = Lim11, lambda = lambda, nu = nu, type = type)) -
          ((Lim21^3) * .E_phiSNI(r = -0.5, a = Lim21, lambda = lambda, nu = nu, type = type) -
             (Lim11^3) * .E_phiSNI(r = -0.5, a = Lim11, lambda = lambda, nu = nu, type = type))
      ) -
      tau * .L(s = 4, lambda = lambda) * (
        2 * (.E_phi(r = -2, a = b_alpha, nu = nu, type = type1) -
               .E_phi(r = -2, a = a_alpha, nu = nu, type = type1)) +
          ((b_alpha^2) * .E_phi(r = -1, a = b_alpha, nu = nu, type = type1) -
             (a_alpha^2) * .E_phi(r = -1, a = a_alpha, nu = nu, type = type1))
      )
  }

  sigma <- sqrt(sigma2)

  EUY1 <- mu + sigma * EUX1

  if (type == "SN" || type == "SCN" || (type == "ST" & nu > 2) || (type == "SSL" & nu > 1)) {
    EUY2 <- mu^2 + 2 * mu * sigma * EUX1 + (sigma^2) * EUX2
  }
  if (type == "SN" || type == "SCN" || (type == "ST" & nu > 3) || (type == "SSL" & nu > 1.5)) {
    EUY3 <- mu^3 + 3 * (mu^2) * sigma * EUX1 + 3 * mu * (sigma^2) * EUX2 + (sigma^3) * EUX3
  }
  if (type == "SN" || type == "SCN" || (type == "ST" & nu > 4) || (type == "SSL" & nu > 2)) {
    EUY4 <- mu^4 + 4 * (mu^3) * sigma * EUX1 + 6 * (mu^2) * (sigma^2) * EUX2 + 4 * mu * (sigma^3) * EUX3 + (sigma^4) * EUX4

    S <- (EUY3 - 3 * EUY1 * EUY2 + 2 * EUY1^3) / ((EUY2 - EUY1^2)^(3/2))
    K <- (EUY4 - 4 * EUY1 * EUY3 + 6 * EUY2 * (EUY1^2) - 3 * (EUY1^4)) / (EUY2 - EUY1^2)^2
    CV <- sqrt((EUY2 - EUY1^2)) / EUY1
  }


  if (type != "ST" & type != "SSL") {
    resp <- c(EUY1 = EUY1,
              EUY2 = EUY2,
              EUY3 = EUY3,
              EUY4 = EUY4,
              S = S,
              K = K,
              CV = CV)
  }
  if (type == "ST") {
    if (nu > 1 & nu <= 2) {
      resp <- c(EUY1 = EUY1)
    }
    if (nu > 2 & nu <= 3) {
      resp <- c(EUY1 = EUY1,
                EUY2 = EUY2)
    }
    if (nu > 3 & nu <= 4) {
      resp <- c(EUY1 = EUY1,
                EUY2 = EUY2,
                EUY3 = EUY3)
    }
    if (nu > 4) {
      resp <- c(EUY1 = EUY1,
                EUY2 = EUY2,
                EUY3 = EUY3,
                EUY4 = EUY4,
                S = S,
                K = K,
                CV = CV)
    }
  }
  if (type == "SSL") {
    if (nu > 0.5 & nu <= 1) {
      resp <- c(EUY1 = EUY1)
    }
    if (nu > 1 & nu <= 1.5) {
      resp <- c(EUY1 = EUY1,
                EUY2 = EUY2)
    }
    if (nu > 1.5 & nu <= 2) {
      resp <- c(EUY1 = EUY1,
                EUY2 = EUY2,
                EUY3 = EUY3)
    }
    if (nu > 2) {
      resp <- c(EUY1 = EUY1,
                EUY2 = EUY2,
                EUY3 = EUY3,
                EUY4 = EUY4,
                S = S,
                K = K,
                CV = CV)
    }
  }

  return(resp)

}

# Estimation ---------------------------
.eSNI <- function(x, mu, sigma2, lambda, nu, a, b, type, shape) {

  pb <- progress::progress_bar$new(
    format = "  processing [:bar]",
    total = 100, clear = TRUE)

  if (type == "SN") {
    LL1 <- function(par) {

      mu <- par[1]
      sigma2 <- par[2]
      lambda <- par[3]

      pdf <- .pdfSNI(y = x, mu = mu, sigma2 = sigma2, lambda = lambda, nu = NULL, type = "SN")/
        (.cdfSNI(y = b, mu = mu, sigma2 = sigma2, lambda = lambda, nu = NULL, type = "SN") -
           .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = NULL, type = "SN"))

      for (i in seq_along(pdf)) {
        if (pdf[i] <= 1e-12) {
          pdf[i] <- 1e-12
        }
      }

      -sum(log(pdf))
    }
    LL1G <- function(par) {

      mu <- par[1]
      sigma2 <- par[2]
      lambda <- par[3]

      try(pb$tick(), silent=TRUE)

      numDeriv::grad(func = LL1,
                     x = c(mu, sigma2, lambda),
                     method.args=list(r = 6))

    }
    fit <- nlminb(start = c(mu = mu, sigma2 = sigma2, lambda = lambda),
                  objective = LL1,
                  gradient = LL1G,
                  lower = c(-1e+12, 0.01, -1e+12),
                  upper = c(1e+12, 1e+12, 1e+12))

  }
  if (type == "ST" || type == "SSL") {
    LL2 <- function(par) {

      mu <- par[1]
      sigma2 <- par[2]
      lambda <- par[3]
      if (shape == FALSE) {
        nu <- par[4]
      }

      pdf <- .pdfSNI(y = x, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type)/
        (.cdfSNI(y = b, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type) -
           .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type))

      for (i in seq_along(pdf)) {
        if (pdf[i] <= 1e-12) {
          pdf[i] <- 1e-12
        }
      }

      -sum(log(pdf))

    }
    LL2G <- function(par) {

      mu <- par[1]
      sigma2 <- par[2]
      lambda <- par[3]
      aux <- c(mu, sigma2, lambda)
      if (shape == FALSE) {
        nu <- par[4]
        aux <- c(mu, sigma2, lambda, nu)
      }

      try(pb$tick(), silent=TRUE)



      numDeriv::grad(func = LL2,
                     x = aux,
                     method.args=list(r = 6))

    }

    if (shape == TRUE) {
      start1 <- c(mu = mu, sigma2 = sigma2, lambda = lambda)
      lower1 <- c(-1e+12, 0.01, -1e+12)
      upper1 <- c(1e+12, 1e+12, 1e+12)
    }
    if (shape == FALSE) {
      start1 <- c(mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu)
      lower1 <- c(-1e+12, 0.01, -1e+12, 1.1)
      upper1 <- c(1e+12, 1e+12, 1e+12, 30)
    }

    fit <- nlminb(start = start1,
                  objective = LL2,
                  gradient = LL2G,
                  lower = lower1,
                  upper = upper1)

  }
  if (type == "SCN") {
    LL3 <- function(par) {

      mu <- par[1]
      sigma2 <- par[2]
      lambda <- par[3]
      nu1 <- ifelse(shape == TRUE, nu[1], par[4])
      nu2 <- ifelse(shape == TRUE, nu[2], par[5])

      pdf <- .pdfSNI(y = x, mu = mu, sigma2 = sigma2, lambda = lambda, nu = c(nu1, nu2), type = type)/
        (.cdfSNI(y = b, mu = mu, sigma2 = sigma2, lambda = lambda, nu = c(nu1, nu2), type = type) -
           .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = c(nu1, nu2), type = type))

      for (i in seq_along(pdf)) {
        if (pdf[i] <= 1e-12) {
          pdf[i] <- 1e-12
        }
      }

      -sum(log(pdf))

    }
    LL3G <- function(par) {

      mu <- par[1]
      sigma2 <- par[2]
      lambda <- par[3]

      aux <- c(mu, sigma2, lambda)
      if (shape == FALSE) {
        nu1 <- par[4]
        nu2 <- par[5]
        aux <- c(mu, sigma2, lambda, nu1, nu2)
      }

      try(pb$tick(), silent=TRUE)

      numDeriv::grad(func = LL3,
                     x = aux,
                     method.args=list(r = 6))

    }

    if (shape == TRUE) {
      start1 <- c(mu = mu, sigma2 = sigma2, lambda = lambda)
      lower1 <- c(-1e+12, 0.01, -1e+12)
      upper1 <- c(1e+12, 1e+12, 1e+12)
    }
    if (shape == FALSE) {
      start1 <- c(mu = mu, sigma2 = sigma2, lambda = lambda, nu1 = nu[1], nu2 = nu[2])
      lower1 <- c(-1e+12, 0.01, -1e+12, 0.01, 0.01)
      upper1 <- c(1e+12, 1e+12, 1e+12, 0.99, 0.99)
    }

    fit <- nlminb(start = start1,
                  objective = LL3,
                  gradient = LL3G,
                  lower = lower1,
                  upper = upper1)

  }

  resp <- fit$par

  cat("\n")
  return(resp)

}

# Generation ---------------------------
.rSNI <- function(n, mu, sigma2, lambda, nu, a, b, type) {

  if (type == "SN") {
    u <- runif(n = n)
    A <- (.cdfSNI(y = b, mu = mu, sigma2 = sigma2, lambda = lambda, nu = NULL, type = type)
          - .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = NULL, type = type))
    aux <- (u * A) + .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type)
    am.x <- sn::qsn(p = aux, xi = mu, omega = sqrt(sigma2), alpha = lambda)
  }
  if (type == "ST") {
    u <- runif(n = n)
    A <- (.cdfSNI(y = b, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type)
          - .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type))
    aux <- (u * A) + .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type)
    am.x <- sn::qst(p = aux, xi = mu, omega = sqrt(sigma2), alpha = lambda, nu = nu)
  }
  if (type == "SSL") {
    u <- rbeta(n = n, shape1 = nu, shape2 = 1)
    a1 <- (a - mu) * sqrt(x = u)
    b1 <- (b - mu) * sqrt(x = u)
    aux0 <- am.x <- c()
    for (i in 1:n) {
      aux0[i] <- .rSNI(n = 1, mu = 0, sigma2 = sigma2, lambda = lambda, nu = NULL, a = a1[i], b = b1[i], type = "SN")
      am.x[i] <- mu + (u[i])^(-1/2) * aux0[i]
    }
  }
  if (type == "SCN") {
    p <- runif(n = n)
    u <- rep(1, n)
    u[p < nu[1]] <- nu[2]
    a1 <- (a - mu) * sqrt(x = u)
    b1 <- (b - mu) * sqrt(x = u)
    aux0 <- am.x <- c()
    for(i in 1:n) {
      aux0[i] <- .rSNI(n = 1, mu = 0, sigma2 = sigma2, lambda = lambda, nu = NULL, a = a1[i], b = b1[i], type = "SN")
      am.x[i] <- mu + (u[i])^(-1/2)*aux0[i]
    }
  }

  return(am.x)

}

# Pdf ---------------------------
.pdfSNI <- function(y, mu, sigma2, lambda, nu, type) {

  z <- (y - mu) / sqrt(sigma2)

  if (type == "SN") {
    resp <- 2 * dnorm(x = y, mean = mu, sd = sqrt(sigma2)) * pnorm(lambda * z)
  }
  if (type == "ST") {
    resp <- 2 * dt(x = z, df = nu) *
      pt(q = sqrt((nu + 1) / (nu + z^2)) * lambda * z, df = nu + 1) / sqrt(sigma2)
  }
  if (type == "SSL") {
    resp <- vector(mode = "numeric", length = length(y))
    for (i in 1:length(y)) {
      f <- function(u) {
        2 * nu * u^(nu - 1) *
          dnorm(x = y[i], mean = mu, sd = sqrt(sigma2/u)) *
          pnorm(q = sqrt(u) * (lambda * z[i]))
      }
      resp[i] <- integrate(Vectorize(f), 0, 1)$value
    }
  }
  if (type == "SCN") {
    resp <- 2 * (nu[1] * dnorm(x = y, mean = mu, sd = sqrt(sigma2/nu[2])) *
                   pnorm(q = sqrt(nu[2]) * lambda * z) +
                   (1 - nu[1]) * dnorm(x = y, mean = mu, sd = sqrt(sigma2)) *
                   pnorm(q = lambda * z))
  }

  return(resp)

}

.dSNI <- function(y, mu, sigma2, lambda, nu, type, a, b) {

  resp <- .pdfSNI(y = y, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type)/
    (.cdfSNI(y = b, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type) -
       .cdfSNI(y = a, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, type = type))

  return(resp)

}

# Cdf ---------------------------
.cdfSNI <- function(y, mu, sigma2, lambda, nu, type) {

  delta <- lambda/sqrt(1 + lambda*lambda)
  Sigma <- matrix(c(1, -delta, -delta, 1), nrow = 2, ncol = 2, byrow = TRUE)
  e1 <- c(1, 0)
  z <- (y-mu)/sqrt(sigma2)

  if (type == "SN") {
    resp <- 2 * mvtnorm::pmvnorm(lower = -Inf, upper = z * e1, sigma = Sigma)[1]
  }
  if (type == "ST") {
    resp <- sn::pst(x = y, xi = mu, omega = sqrt(sigma2), alpha = lambda, nu = nu)
  }
  if (type == "SSL") {
    f <- function(u) {
      resp <-
        2 * nu *
        mvtnorm::pmvnorm(lower = -Inf, upper = sqrt(u) * y *e1, sigma = Sigma)[1] *
        u ^ (nu - 1)
    }
    resp <- integrate(Vectorize(f), 0, 1)$value
  }
  if (type == "SCN") {
    resp <-
      2 * (nu[1] * mvtnorm::pmvnorm(lower = -Inf, upper = sqrt(nu[2]) * y * e1, sigma = Sigma)[1] +
             (1 - nu[1]) * mvtnorm::pmvnorm(lower = -Inf, upper = y * e1, sigma = Sigma)[1])
  }

  return(resp)

}

# Extras ---------------------------
.GamaInc <- function(a, b) {

  f <- function(t) {
    exp(-t) * t^(a-1)
  }

  resp <- integrate(f, 0, b)$value

  return(resp)

}

.L <- function(s, lambda) {
  (sqrt(2/pi) * lambda) / ((1 + lambda^2)^(s/2))
}

.Bounds <- function(mu, sigma2, lambda, nu, a, b, type) {

  delta <- lambda / sqrt(1 + lambda^2)
  Delta <- sqrt(sigma2) * delta

  if (lambda == 0) {
    if (mu < a || mu > b) {
      stop("mu must be between lower and upper bounds.")
    }
  } else {
    if (type == "SN") {
      lim <- sqrt(2/pi) * Delta
      if (mu < a - lim|| mu > b - lim) {
        stop("check the upper and lower bounds for asymmetric cases.")
      }
    }
    if (type == "ST") {
      lim <- sqrt(2/pi) * sqrt(nu/2) * (gamma(0.5*(nu - 1))/gamma(0.5*nu)) * Delta
      if (mu < a - lim|| mu > b - lim) {
        stop("check the upper and lower bounds for asymmetric cases.")
      }
    }
    if (type == "SSL") {
      lim <- sqrt(2/pi) * (2*nu/(2*nu - 1))  * Delta
      if (mu < a - lim|| mu > b - lim) {
        stop("check the upper and lower bounds for asymmetric cases.")
      }
    }
    if (type == "SCN") {
      lim <- sqrt(2/pi) * ((nu[1]/sqrt(nu[2])) + 1 - nu[1])  * Delta
      if (mu < a - lim|| mu > b - lim) {
        stop("check the upper and lower bounds for asymmetric cases.")
      }
    }
  }

}

# E_phi and E_Phi ---------------------------
.E_phi <- function(r, a, nu, delta, type) {

  if (type == "N") {
    resp <- dnorm(a)
  }
  if (type == "T") {
    Aux0 <- gamma((nu + 2*r) / 2)
    Aux1 <- gamma(nu / 2) * sqrt(2 * pi)
    Aux2 <- Aux0 / Aux1
    Aux3 <- (nu / 2) ^ (nu / 2)
    Aux4 <- (((a^2) + nu) / 2) ^ (-(nu + 2*r)/2)
    resp <- Aux2*Aux3*Aux4
  }
  if (type == "SL") {
    Aux0 <- nu/sqrt(x = 2*pi)
    Aux1 <- (0.5*(a^2))^(-(nu+r))
    Aux2 <- .GamaInc(a = nu + r, b = 0.5*(a^2))
    resp <- Aux0*Aux1*Aux2
  }
  if (type == "CN") {
    Aux0 <- nu[1] * (nu[2]^(r)) * dnorm(x = a*sqrt(nu[2]))
    Aux1 <- (1 - nu[1]) * dnorm(x = a)
    resp <- Aux0 + Aux1
  }

  return(resp)

}

.E_Phi <- function(r, a, nu, delta, type) {

  if (type == "N") {
    resp <- pnorm(a)
  }
  if (type == "T") {
    Aux0 <- gamma(0.5*nu + r)
    Aux1 <- gamma(nu/2)
    Aux2 <- Aux0/Aux1
    Aux3 <- (2/nu)^(r)
    Aux4 <- sqrt((2*r + nu)/nu) * a
    Aux5 <- pt(Aux4, df = 2*r + nu)
    resp <- Aux2*Aux3*Aux5
  }
  if (type == "SL") {
    Aux0 <- nu/(nu+r)
    Aux1 <- .cdfSNI(y = a, mu = 0, sigma2 = 1, lambda = 0, nu = nu + r, type = "SSL")
    resp <- Aux0*Aux1
  }
  if (type == "CN") {
    Aux0 <- nu[2]^(r) *
      .cdfSNI(y = a, mu = 0, sigma2 = 1, lambda = 0, nu = nu, type = "SCN")
    Aux1 <- (1 - nu[1]) * (1 - nu[2]^r) * pnorm(a)
    resp <- Aux0 + Aux1
  }

  return(resp)

}

# E_phiSNI and E_PhiSNI ---------------------------
.E_phiSNI <- function(r, a, nu, delta, lambda, type) {

  if (type == "SN") {
    resp <- .pdfSNI(y = a, mu = 0, sigma2 = 1, lambda = lambda, nu  = NULL, type = "SN")
  }
  if (type == "ST") {
    Aux0 <- (2^(r+1)) * (nu^(nu/2)) * gamma((nu+2*r)/2)
    Aux1 <- sqrt(2*pi) * gamma(nu/2) * ((a^2)+nu)^((nu+2*r)/2)
    Aux2 <- Aux0/Aux1
    Aux3 <- sqrt((2*r+nu)/(a^2+nu))*lambda*a
    Aux5 <- pt(Aux3, df = (2*r+nu))
    resp <- Aux2*Aux5
  }
  if (type == "SSL") {
    f <- function(u) {
      .pdfSNI(y = a*sqrt(u), mu = 0, sigma2 = 1, lambda = lambda, nu = NULL, type = "SN") *
        nu * (u^(nu - 1)) *
        (u^r)
    }
    resp <- integrate(f, 0, 1)$value
  }
  if (type == "SCN") {
    Aux0 <- nu[1] * (nu[2]^(r)) *
      .pdfSNI(y = a * sqrt(nu[2]), mu = 0, sigma2 = 1, lambda = lambda, nu = NULL, type = "SN")
    Aux1 <- (1 - nu[1]) *
      .pdfSNI(y = a, mu = 0, sigma2 = 1, lambda = lambda, nu = NULL, type = "SN")
    resp <- Aux0 + Aux1
  }

  return(resp)

}

.E_PhiSNI <- function(r, a, nu, delta, lambda, type) {

  if (type == "SN") {
    resp <- .cdfSNI(y = a, mu = 0, sigma2 = 1, lambda = lambda, nu = NULL, type = "SN")
  }
  if (type == "ST") {

    Aux0 <- gamma(0.5*nu + r)/gamma(0.5 * nu)
    Aux1 <- (2/nu) ^ r

    alpha1 <- 0.5*nu + r
    alpha2 <- 0.5*nu

    f <- function(u) {
      .cdfSNI(y = a*sqrt(u), mu = 0, sigma2 = 1, lambda = lambda, nu = NULL, type = "SN") *
        (1/gamma(alpha1)) * (alpha2 ^ alpha1) * (u ^ (alpha1 - 1)) * exp(-u * alpha2)
    }

    Aux2 <- integrate(Vectorize(f), 0, Inf)$value

    resp <- Aux0 * Aux1 * Aux2

  }
  if (type == "SSL") {
    f <- function(u) {
      (u^r) *
        .cdfSNI(y = a*sqrt(u), mu = 0, sigma2 = 1, lambda = lambda, nu = NULL, type = "SN") *
        nu * (u^(nu - 1))
    }
    resp <- integrate(Vectorize(f), 0, 1)$value
  }
  if (type == "SCN") {

    delta <- lambda/sqrt(1 + lambda*lambda)
    Sigma <- matrix(c(1, -delta, -delta, 1), nrow = 2, ncol = 2, byrow = TRUE)
    e1 <- c(1, 0)

    Aux0 <- nu[2]^(r) *
      .cdfSNI(y = a, mu = 0, sigma2 = 1, lambda = lambda, nu = nu, type = "SCN")
    Aux1 <- 2 * (1 - nu[1]) * (1 - nu[2]^r) *
      mvtnorm::pmvnorm(lower = -Inf, upper = a * e1, sigma = Sigma)
    resp <- Aux0 + Aux1
  }

  return(resp)

}

