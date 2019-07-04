#' Estimate the parameters of Truncated Scale Mixtures of Skew-Normal Distributions
#'
#' This function obtains the maximum likelihood estimators of the TSMSN (Skew-Normal, Skew-t, Skew-Slash and Skew-Contaminated Normal) distribution parameters by direct maximization.
#'
#' @param x Dataset.
#' @param mu Initial location parameter (Optional).
#' @param sigma2 Initial scale parameter (Optional).
#' @param lambda Initial skewness parameter (Optional).
#' @param nu Fixed shape parameter.  Must be NULL in case of Skew-Normal distribution.  Must be a bidimensional vector in case of skew-contaminated normal distribution (SCN) and contaminated normal distribution (CN).
#' @param a Lower bound.
#' @param b Upper bound.
#' @param dist Distribution to be used:  "SN" for Skew-Normal model, "ST" for Skew-t model, "SSL" for Skew-slash model and "SCN" for Skew-contaminated Normal model.
#' @param shape For ST, SSL and SCN distribution. Consider the parameter nu as fixed and known. If TRUE nu must be provided.
#'
#' @details For the SMN family, consider lambda = 0. For the Skew-contaminated Normal and Contaminated Normal distribution, each component of the bidimensional vector "nu" must be on (0,1). For the estimation in the cases of distributions ST, SSL and SCN nu is considered fixed, but may be known or unknown. The shape parameter is the one that regulates if nu is known or unknown.
#'
#' @return Returns a vector with the maximum likelihood estimators of the distribution parameters.
#'
#' @references
#'
#' Lachos, V. H.; Garay, A. M.; Cabral, C. R. "Moments of truncated scale mixtures of skew-normal distributions." Brazilian Journal of Probability and Statistics (In press).
#'
#' Basso, Rodrigo M., et al. "Robust mixture modeling based on scale mixtures of skew-normal distributions." Computational Statistics & Data Analysis 54.12 (2010): 2926-2941.
#'
#' @importFrom stats dnorm dt integrate pnorm pt rbeta runif optim nlminb var
#'
#' @export
#'
#' @examples
#' x <- rTSMSN(n = 100, mu = 0, sigma2 = 1, lambda = 0, nu = NULL, a = -2, b = 2, dist = "SN")
#' eTSMSN(x, a = -2, b = 2, dist = "SN")
eTSMSN <- function(x, mu = 0.01, sigma2 = 1.01, lambda = 0, nu = 5, a = -Inf, b = Inf, dist = "SN", shape = FALSE) {

  if (!is.logical(shape)) {
    stop("shape must be logical.")
  }

  if (!is.null(nu)) {
    if (length(nu) == 1) {
      if (nu == 5 & dist == "SCN") {
        nu <- c(0.5, 0.5)
      }
    }
  }

  if (length(x) < 1)
    stop("x parameter must be a vector.")

  if (length(mu) > 1)
    stop("mu parameter must be a scalar.")

  if (length(sigma2) > 1)
    stop("sigma2 parameter must be a scalar.")
  if (sigma2 <= 0)
    stop("Scale parameter must be positive.")

  if (length(lambda) > 1)
    stop("lambda parameter must be a scalar.")

  if (a > b)
    stop("a must be lower than b.")
  if (length(a) > 1 || length(b) > 1)
    stop("Range values must be scalar.")

  if ((dist != "SN") && (dist != "ST") && (dist != "SSL") && (dist != "SCN"))
    stop("Distribution family not supported. Check documentation.")

  if (dist == "SCN") {
    if (length(nu) != 2)
      stop("nu must be a bidimensional vector in case of Contaminated Normal distribution.")
    if (nu[1] <= 0 || nu[1] >= 1)
      stop("nu[1] must lies in (0,1).")
    if (nu[2] <= 0 || nu[2] >= 1)
      stop("nu[2] must lies in (0,1).")
  }

  if (a < -1e+12) {
    a <- -1e+12
  }
  if (b > 1e+12) {
    b <- 1e+12
  }

  if (mu > 1e+12) {
    mu <- 1e+12
  }
  if (mu < -1e+12) {
    mu <- -1e+12
  }
  if (sigma2 > 1e+12) {
    sigma2 <- 1e+12
  }
  if (sigma2 < -1e+12) {
    sigma2 <- -1e+12
  }
  if (lambda > 1e+12) {
    lambda <- 1e+12
  }
  if (lambda < -1e+12) {
    lambda <- -1e+12
  }

  if (dist != "SCN" && dist != "SN"){
    if (length(nu) > 1)
      stop("nu parameter must be a scalar.")
    if (length(nu) == 0)
      stop("initial value for nu parameter must be provided in case of ST and SSL.")
    if (dist == "ST") {
      if (nu <= 1)
        stop("nu must be greater than 1 for ST.")
    }
    if (dist == "SSL") {
      if (nu <= 0.5)
        stop("nu must be greater than 0.5 for SSL.")
    }
    if (nu >= 30) {
      nu = 30
    }
  }

  if (mu == 0.01) {
    mu <- mean(x, na.rm = TRUE)
  }
  if (sigma2 == 1.01) {
    sigma2 <- var(x, na.rm = TRUE)
  }

  .Bounds(mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, a = a, b = b, type = dist)

  .eSNI(x = x, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, a = a, b = b, type = dist, shape = shape)

}

#' Generate the Truncated Scale Mixtures of Skew-Normal Distributions
#'
#' This function generate random samples from Skew-Normal, Skew-t, Skew-Slash and Skew-Contaminated Normal, using the the inverse method.
#'
#' @param n Number of observations.
#' @param mu Location parameter.
#' @param sigma2 Scale parameter.
#' @param lambda Skewness parameter.
#' @param nu Shape parameter.  Must be NULL in case of Skew-Normal distribution.  Must be a bidimensional vector in case of skew-contaminated normal distribution (SCN) and contaminated normal distribution (CN).
#' @param a Lower bound.
#' @param b Upper bound.
#' @param dist Distribution to be used:  "SN" for Skew-Normal model, "ST" for Skew-t model, "SSL" for Skew-slash model and "SCN" for Skew-contaminated Normal model.
#'
#' @details For the SMN family, consider lambda = 0. For the Skew-contaminated Normal and Contaminated Normal distribution, each component of the bidimensional vector "nu" must be on (0,1).
#'
#' @return Returns a vector with the sample generated according to the distribution and parameters defined.
#'
#' @references
#'
#' Lachos, V. H.; Garay, A. M.; Cabral, C. R. "Moments of truncated scale mixtures of skew-normal distributions." Brazilian Journal of Probability and Statistics (In press).
#'
#' Basso, Rodrigo M., et al. "Robust mixture modeling based on scale mixtures of skew-normal distributions." Computational Statistics & Data Analysis 54.12 (2010): 2926-2941.
#'
#' @importFrom stats dnorm dt integrate pnorm pt rbeta runif optim nlminb var
#'
#' @export
#'
#' @examples
#' rTSMSN(n = 100, mu = 0, sigma2 = 1, nu = NULL, lambda = 0, a = -Inf, b = Inf, dist = "SN")
rTSMSN <- function(n, mu = 0, sigma2 = 1, lambda = 0, nu = NULL, a = -Inf, b = Inf, dist = "SN") {

  if (length(n) > 1)
    stop("n parameter must be a scalar.")
  if (n <= 0)
    stop("n parameter must be positive.")

  if (length(mu) > 1)
    stop("mu parameter must be a scalar.")

  if (length(sigma2) > 1)
    stop("sigma2 parameter must be a scalar.")
  if (sigma2 <= 0)
    stop("Scale parameter must be positive.")

  if (length(lambda) > 1)
    stop("lambda parameter must be a scalar.")

  if (a > b)
    stop("a must be lower than b.")
  if (length(a) > 1 || length(b) > 1)
    stop("Range values must be scalar.")

  if ((dist != "SN") && (dist != "ST") && (dist != "SSL") && (dist != "SCN"))
    stop("Distribution family not supported. Check documentation.")

  if (dist == "SCN") {
    if (length(nu) != 2)
      stop("nu must be a bidimensional vector in case of Contaminated Normal distribution.")
    if (nu[1] <= 0 || nu[1] >= 1)
      stop("nu[1] must lies in (0,1).")
    if (nu[2] <= 0 || nu[2] >= 1)
      stop("nu[2] must lies in (0,1).")
  }

  if (a < -1e+12) {
    a <- -1e+12
  }
  if (b > 1e+12) {
    b <- 1e+12
  }

  if (mu > 1e+12) {
    mu <- 1e+12
  }
  if (mu < -1e+12) {
    mu <- -1e+12
  }
  if (sigma2 > 1e+12) {
    sigma2 <- 1e+12
  }
  if (sigma2 < -1e+12) {
    sigma2 <- -1e+12
  }
  if (lambda > 1e+12) {
    lambda <- 1e+12
  }
  if (lambda < -1e+12) {
    lambda <- -1e+12
  }

  if (dist != "SCN" && dist != "SN"){
    if (length(nu) > 1)
      stop("nu parameter must be a scalar.")
    if (length(nu) == 0)
      stop("initial value for nu parameter must be provided in case of ST and SSL.")
    if (dist == "ST") {
      if (nu <= 1)
        stop("nu must be greater than 1 for ST.")
    }
    if (dist == "SSL") {
      if (nu <= 0.5)
        stop("nu must be greater than 0.5 for SSL.")
    }
    if (nu >= 30) {
      nu = 30
    }
  }

  .Bounds(mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, a = a, b = b, type = dist)

  .rSNI(n = n, mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, a = a, b = b, type = dist)

}

#' Moments of Truncated Scale Mixtures of Skew-Normal Distributions
#'
#' Return the first four moments of the TSMSN distributions (Skew Normal, Skew t, Skew Slash or Skew Contaminated Normal).
#'
#' @param mu Location parameter.
#' @param sigma2 Scale parameter.
#' @param lambda Skewness parameter.
#' @param nu Shape parameter.  Must be NULL in case of Skew-Normal distribution.  Must be a bidimensional vector in case of skew-contaminated normal distribution (SCN) and contaminated normal distribution (CN).
#' @param a Lower bound.
#' @param b Upper bound.
#' @param dist Distribution to be used:  "SN" for Skew-Normal model, "ST" for Skew-t model, "SSL" for Skew-slash model and "SCN" for Skew-contaminated Normal model.
#' @param empir If TRUE provides the empirical moments.
#'
#' @details For the SMN family, consider lambda = 0. For the Skew-contaminated Normal and Contaminated Normal distribution, each component of the bidimensional vector "nu" must be on (0,1).
#'
#' @return Returns the four moments, the skewness (S), kurtosis (k) and coefficient of variation(CV). If “empir = TRUE”, returns also the Empirical moments.
#'
#' @references
#'
#' Lachos, V. H.; Garay, A. M.; Cabral, C. R. "Moments of truncated scale mixtures of skew-normal distributions." Brazilian Journal of Probability and Statistics (In press).
#'
#' Basso, Rodrigo M., et al. "Robust mixture modeling based on scale mixtures of skew-normal distributions." Computational Statistics & Data Analysis 54.12 (2010): 2926-2941.
#'
#' @importFrom stats dnorm dt integrate pnorm pt rbeta runif optim nlminb var
#'
#' @export
#'
#' @examples
#' mTSMSN(mu = 1, sigma2 = 1, nu = NULL, lambda = 1, a = -2, b = 2, dist = "SN", empir = TRUE)
mTSMSN <- function(mu = 0, sigma2 = 1, lambda = 0, nu = NULL, a = -Inf, b = Inf, dist = "SN", empir = TRUE) {

  if (length(mu) > 1)
    stop("mu parameter must be a scalar.")

  if (length(sigma2) > 1)
    stop("sigma2 parameter must be a scalar.")
  if (sigma2 <= 0)
    stop("Scale parameter must be positive.")

  if (length(lambda) > 1)
    stop("lambda parameter must be a scalar.")

  if (a > b)
    stop("a must be lower than b.")
  if (length(a) > 1 || length(b) > 1)
    stop("Range values must be scalar.")

  if ((dist != "SN") && (dist != "ST") && (dist != "SSL") && (dist != "SCN"))
    stop("Distribution family not supported. Check documentation.")

  if (dist == "SCN") {
    if (length(nu) != 2)
      stop("nu must be a bidimensional vector in case of Contaminated Normal distribution.")
    if (nu[1] <= 0 || nu[1] >= 1)
      stop("nu[1] must lies in (0,1).")
    if (nu[2] <= 0 || nu[2] >= 1)
      stop("nu[2] must lies in (0,1).")
  }

  if (a < -1e+12) {
    a <- -1e+12
  }
  if (b > 1e+12) {
    b <- 1e+12
  }

  if (mu > 1e+12) {
    mu <- 1e+12
  }
  if (mu < -1e+12) {
    mu <- -1e+12
  }
  if (sigma2 > 1e+12) {
    sigma2 <- 1e+12
  }
  if (sigma2 < -1e+12) {
    sigma2 <- -1e+12
  }
  if (lambda > 1e+12) {
    lambda <- 1e+12
  }
  if (lambda < -1e+12) {
    lambda <- -1e+12
  }

  if (dist != "SCN" && dist != "SN"){
    if (length(nu) > 1)
      stop("nu parameter must be a scalar.")
    if (length(nu) == 0)
      stop("initial value for nu parameter must be provided in case of ST and SSL.")
    if (dist == "ST") {
      if (nu <= 1)
        stop("nu must be greater than 1 for ST.")
    }
    if (dist == "SSL") {
      if (nu <= 0.5)
        stop("nu must be greater than 0.5 for SSL.")
    }
    if (nu >= 30) {
      nu = 30
    }
  }

  if (!is.logical(empir)) {
    stop("empir must be logical.")
  }

  .Bounds(mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, a = a, b = b, type = dist)

  if (empir == TRUE) {
    x <- TSMSN::rTSMSN(n = 20000, mu = mu, sigma2 = sigma2, nu = nu, lambda = lambda, a = a, b = b, dist = dist)
    resp_empir <- round(c(EUY1 = mean(x), EUY2 = mean(x^2), EUY3 = mean(x^3), EUY4 = mean(x^4)),
                        digits = 5)
  }

  resp_theor <- round(.mSNI(mu = mu, sigma2 = sigma2, lambda = lambda, nu = nu, a = a, b = b, type = dist),
                digits = 5)

  if (length(resp_theor) != 7) {
    resp_empir <- resp_empir[1:length(resp_theor)]
  }

  if (empir == TRUE) {
    resp <- list("Theoretical Moments" = resp_theor,
                 "Empirical Moments" = resp_empir)
  }
  if (empir == FALSE) {
    resp <- resp_theor
  }

  return(resp)

}
