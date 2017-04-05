#' Function for estimating tuning parameters
#'
#' \code{estRegPars}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @param \code{delta} ridge regression parameter for when X is not full rank
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat} and \code{kappa.hat}
#' @export
estRegPars <-function(y, X, delta = 0) {


  p <- ncol(X)
  n <- nrow(X)

  XtX <- crossprod(X)

  full.X <- min(eigen(XtX)$values) > 0
  if (!full & delta == 0) {delta <- 1}

  D <- solve(XtX + delta^2*diag(p))
  A <- diag(n) - tcrossprod(tcrossprod(X, D), X)
  b <- crossprod(D, crossprod(X, y))
  r <- crossprod(A, y)

  DXtX <- crossprod(D, XtX)
  XtXDDXtX <- crossprod(DXtX)
  DXt <- tcrossprod(D, X)
  XDDXt <- crossprod(DXt)
  AX <- crossprod(A, X)

  e.bb <- sum(diag(XtXDDXtX))/p
  e.be <- sum(diag(XDDXt))/p
  e.eb <- sum(diag(crossprod(AX)))/n
  e.ee <- sum(diag(crossprod(A)))/n

  E <- rbind(c(e.bb, e.be),
             c(e.eb, e.ee))

  sig.2.ests <- solve(E)%*%c(mean(b^2), mean(r^2))
  sigma.beta.sq.hat <- ifelse(sig.2.ests[1] > 0, sig.2.ests[1], 0)
  sigma.epsi.sq.hat <- ifelse(sig.2.ests[2] > 0, sig.2.ests[2], 0)

  aa <- sum(3*(diag(XtXDDXtX)*sigma.beta.sq.hat + diag(XDDXt)*sigma.epsi.sq.hat)^2/p)
  bb <- sum(diag(DXtX^4))/p
  kappa.hat <- (mean(b^4) - aa)/(bb*sigma.beta.sq.hat^2)

  return(list = c("sigma.beta.sq.hat" = sigma.beta.sq.hat,
                  "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
                  "kappa.hat" = kappa.hat))

}

