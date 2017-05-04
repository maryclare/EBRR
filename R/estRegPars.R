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

  aa <- sum(3*(diag(XtHXDDXtHX)*sigma.beta.sq.hat + diag(DXtHXD)*sigma.epsi.sq.hat)^2/p)
  bb <- sum(diag(DXtHX^4))/p
  kappa.hat <- (mean(b^4) - aa)/(bb*sigma.beta.sq.hat^2)

  return(list = c("sigma.beta.sq.hat" = sigma.beta.sq.hat,
                  "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
                  "kappa.hat" = kappa.hat))

}


