fq <- function(q, kurt) {
  gamma(5/q)*gamma(1/q)/(gamma(3/q)^2) - kurt
}

fpq <- function(q) {
  (gamma(5/q)*gamma(1/q)/(q^2*gamma(3/q)^2))*(6*digamma(3/q) - digamma(1/q) - 5*digamma(5/q))
}

# Use Newton's method: https://en.wikipedia.org/wiki/Newton%27s_method
nrq <- function(kurt, sval = 0.032, tol = 10^(-12)) { # This starting value is the lowest possible
  if (kurt < 6) {
    sval <- 1
  } else if (kurt < 3) {
    sval <- 2
  }
  x.old <- Inf
  x.new <- sval
  while (abs(x.new - x.old) > tol) {
    x.old <- x.new
    x.new <- x.old - fq(x.old, kurt)/fpq(x.old)
  }
  return(x.new)
}
#' Function for estimating tuning parameters
#'
#' \code{estRegPars}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @param \code{delta} ridge regression parameter for when X is not full rank
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat} and \code{kappa.hat}
#' @export
estRegPars <-function(y, X, W = NULL, delta = 0, precomp = NULL) {


  p <- ncol(X)
  n <- nrow(X)

  if (is.null(W) & is.null(precomp)) {
    H <- diag(n)
  } else if (!is.null(W) & is.null(precomp)) {
    H <- diag(n) - tcrossprod(tcrossprod(W, solve(crossprod(W))),W)
  } else if (!is.null(precomp)) {
    H <- precomp[["H"]]
  }

  if (is.null(precomp)) {
    HX <- crossprod(H, X)
    XtHX <- crossprod(HX)
    full.X <- min(eigen(XtHX)$values) > 0
    if (!full.X & delta == 0) {delta <- 1}
    D <- solve(XtHX + delta^2*diag(p))
    A <- diag(n) - tcrossprod(tcrossprod(HX, D), HX)
    DXtHX <- crossprod(D, XtHX)
    XtHXDDXtHX <- crossprod(DXtHX)
    HXD <- tcrossprod(HX, D)
    DXtHXD <- crossprod(HXD)
    AHX <- crossprod(A, HX)
    HA <- crossprod(H, A)
    XtHAAHX <- crossprod(AHX)
    AHA <- crossprod(HA)
  } else {
    HX <- precomp[["HX"]]
    D <- precomp[["D"]]
    A <- precomp[["A"]]
    DXtHX <- precomp[["DXtHX"]]
    XtHXDDXtHX <- precomp[["XtHXDDXtHX"]]
    DXtHXD <- precomp[["DXtHXD"]]
    XtHAAHX <- precomp[["XtHAAHX"]]
    AHA <- precomp[["AHA"]]
  }

  b <- crossprod(D, crossprod(HX, y))
  r <- crossprod(A, crossprod(H, y))

  e.bb <- sum(diag(XtHXDDXtHX))/p
  e.be <- sum(diag(DXtHXD))/p
  e.eb <- sum(diag(XtHAAHX))/n
  e.ee <- sum(diag(AHA))/n

  E <- rbind(c(e.bb, e.be),
             c(e.eb, e.ee))

  sig.2.ests <- solve(E)%*%c(mean(b^2), mean(r^2))
  sigma.beta.sq.hat <- sig.2.ests[1]
  sigma.epsi.sq.hat <- sig.2.ests[2]

  alpha.beta <- sum(diag(crossprod(DXtHX)))/p
  gamma.beta <- sum(diag(DXtHX^4))/p
  omega.beta <- 3*(sum(diag(crossprod(DXtHX)^2)) - sum(diag(DXtHX^4)))/p

  b <- crossprod(HXD, y)

  test.stat <- (mean(b^4))/(mean(b^2)^2)

  kappa.hat <- (alpha.beta^2/gamma.beta)*(test.stat - omega.beta/alpha.beta^2)
  q.hat <- nrq(kappa.hat)

  # Old "unbiased"
  # aa <- sum(3*(diag(XtHXDDXtHX)*sigma.beta.sq.hat + diag(DXtHXD)*sigma.epsi.sq.hat)^2/p)
  # bb <- sum(diag(DXtHX^4))/p
  # kappa.hat <- (mean(b^4) - aa)/(bb*sigma.beta.sq.hat^2)

  return(list = c("sigma.beta.sq.hat" = sigma.beta.sq.hat,
                  "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
                  "kappa.hat" = kappa.hat,
                  "q.hat" = q.hat,
                  "test.stat" = test.stat))

}


