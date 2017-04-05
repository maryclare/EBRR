#' Function for estimating tuning parameters
#'
#' \code{estRegPars}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @param \code{delta} ridge regression parameter for when X is not full rank
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat} and \code{kappa.hat}
#' @export
estRegPars <-function(y, X, delta = 0, XtX = NULL, D = NULL, A = NULL, DXtX = NULL,
                      XtXDDXtX = NULL, XD = NULL, DXtXD = NULL, AX = NULL) {


  p <- ncol(X)
  n <- nrow(X)

  if (is.null(XtX)) {
    XtX <- crossprod(X)
  }

  full.X <- min(eigen(XtX)$values) > 0
  if (!full.X & delta == 0) {delta <- 1}

  if (is.null(D)) {
    D <- solve(XtX + delta^2*diag(p))
  }
  if (is.null(A)) {
    A <- diag(n) - tcrossprod(tcrossprod(X, D), X)
  }
  b <- crossprod(D, crossprod(X, y))
  r <- crossprod(A, y)

  if (is.null(DXtX)) {
    DXtX <- crossprod(D, XtX)
  }
  if (is.null(XtXDDXtX)) {
    XtXDDXtX <- crossprod(DXtX)
  }
  if (is.null(XD)) {
    XD <- tcrossprod(X, D)
  }
  if (is.null(DXtXD)) {
    DXtXD <- crossprod(XD)
  }
  if (is.null(AX)) {
    AX <- crossprod(A, X)
  }
  e.bb <- sum(diag(XtXDDXtX))/p
  e.be <- sum(diag(DXtXD))/p
  e.eb <- sum(diag(crossprod(AX)))/n
  e.ee <- sum(diag(crossprod(A)))/n

  E <- rbind(c(e.bb, e.be),
             c(e.eb, e.ee))

  sig.2.ests <- solve(E)%*%c(mean(b^2), mean(r^2))
  sigma.beta.sq.hat <- ifelse(sig.2.ests[1] > 0, sig.2.ests[1], 0)
  sigma.epsi.sq.hat <- ifelse(sig.2.ests[2] > 0, sig.2.ests[2], 0)

  aa <- sum(3*(diag(XtXDDXtX)*sigma.beta.sq.hat + diag(DXtXD)*sigma.epsi.sq.hat)^2/p)
  bb <- sum(diag(DXtX^4))/p
  kappa.hat <- (mean(b^4) - aa)/(bb*sigma.beta.sq.hat^2)
  kappa.hat.sigma.beta.qd <- (((p + 2)/(p - 1))*(mean(b^4) - 3*p*(mean(b^2)^2)/(p + 2)))

  return(list = c("sigma.beta.sq.hat" = sigma.beta.sq.hat,
                  "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
                  "kappa.hat" = kappa.hat,
                  "kappa.hat.sigma.beta.qd" = kappa.hat.sigma.beta.qd))

}


