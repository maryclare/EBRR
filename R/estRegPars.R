# Function for estimating variance parameters from likelihood
rrmmle<-function(y,X,emu=FALSE,s20=1,t20=1)
{
  sX<-svd(X, nu = nrow(X), nv = ncol(X))
  lX<-sX$d^2
  tUX<-t(sX$u)
  xs<-apply(X,1,sum)

  if(nrow(X)>ncol(X))
  {
    lX<-c(lX,rep(0,nrow(X)-ncol(X)))
  }

  objective<-function(ls2t2,emu)
  {
    s2<-exp(ls2t2[1]) ; t2<-exp(ls2t2[2])
    mu<-emu*sum((tUX%*%xs)*(tUX%*%y)/(lX*t2+s2))/sum((tUX%*%xs)^2/(lX*t2+s2))
    sum(log( lX*t2 + s2 )) + sum( (tUX%*%(y-mu*xs))^2/(lX*t2+s2) ) # + 1/s2^10 # Keep s2 = 0 from being a mode
  }

  fit<-optim(log(c(s20,t20)),objective,emu=emu, method = "L-BFGS-B", lower = c(-10, -Inf))
  s2<-exp(fit$par[1]) ; t2<-exp(fit$par[2])
  mu<-emu*sum((tUX%*%xs)*(tUX%*%y)/(lX*t2+s2))/sum((tUX%*%xs)^2/(lX*t2+s2))
  if (fit$convergence == 0) {
    return(c(mu,t2,s2))
  } else {
    return(rep(NA, 3))
  }
}
#
# lik <- function(pars, XXt.val, Uty) {
#   rho <- exp(pars[1])
#   sig.sq <- exp(pars[2])
#   v <- (XXt.val*rho + 1)
#
#   as.numeric(sig.sq*(length(v)*log(sig.sq) + sum(log(v))) + sum(Uty^2/(v))) + 1/sig.sq
# }
#
# XXt.val <- svd(X, nu = nrow(X))$d^2
# Uty <- t(svd(X, nu = nrow(X))$u)%*%y
#
# lik.rho <- function(pars, XXt.val, Uty, sig.sq) {
#   rho <- exp(pars[1])
#   v <- (XXt.val*rho + 1)
#
#   as.numeric(sig.sq*(length(v)*log(sig.sq) + sum(log(v))) + sum(Uty^2/(v)))
# }
#
# up.sig <- function(XXt.val, Uty, rho) {
#   mean(Uty^2/(XXt.val*rho + 1))
# }
#
# alternating <- function(XXt.val, Uty, sig.sq.start = 10, iter = 100) {
#   pars <- matrix(nrow = iter + 1, ncol = 3)
#   pars[1, 1] <- sig.sq.start
#   for (i in 2:(iter + 1)) {
#   opt <- optim(par = rnorm(1),
#                fn = lik.rho, XXt.val = XXt.val, Uty = Uty, sig.sq = pars[i - 1, 1], method = "BFGS")
#   pars[i, 2] <- ifelse(opt$convergence == 0, exp(opt$par[1]), NA)
#   pars[i, 1] <- up.sig(XXt.val, Uty, pars[i, 2])
#   }
#   pars <- pars[-1, ]
#
# }
#
#
# sig.sqs <- seq(0.0001, 0.0025, length.out = 5)
# rhos <- seq(0.5, 10, length.out = 1000)
# vals <- matrix(nrow = length(rhos), ncol = length(sig.sqs))
# for (i in 1:length(rhos)) {
#   for (j in 1:length(sig.sqs)) {
#     v <- sig.sqs[j]*(XXt.val*rhos[i] + 1)
#     vals[i, j] <- as.numeric(sum(log(v)) + sum(Uty^2/v))
#   }
# }
# plot(rhos, vals[, 1], ylim = range(vals), type = "l")
# for (j in 2:length(sig.sqs)) {
#   lines(rhos, vals[, j], lty = j)
# }
# legend("topright", lty = 1:length(sig.sqs), legend = sig.sqs)
#
# sig.sqs <- seq(10^(-37), 10^(-35), length.out = 1000)
# rhos <- seq(6.940486e+35, 6.940486e+35, length.out = 1)
# vals <- matrix(nrow = length(rhos), ncol = length(sig.sqs))
# for (i in 1:length(rhos)) {
#   for (j in 1:length(sig.sqs)) {
#     v <- sig.sqs[j]*(XXt.val*rhos[i] + 1)
#     vals[i, j] <- as.numeric(sum(log(v)) + sum(Uty^2/v))
#   }
# }
# plot(sig.sqs, vals[1, ], ylim = range(vals), type = "l")
# for (j in 2:length(rhos)) {
#   lines(sig.sqs, vals[j, ], lty = j)
# }
# legend("topright", lty = 1:length(rhos), legend = rhos)
# lik.grad <- function(pars, XXt.val, Uty) {
#   tau.sq <- exp(pars[1])
#   sig.sq <- exp(pars[2])
#   v <- XXt.val*tau.sq + sig.sq
#
#   c(as.numeric(sum(XXt.val/v) - sum(Uty^2*XXt.val/v^2)),
#     as.numeric(sum(1/v) - sum(Uty^2/v^2)))
# }
#
#
# exp(optim(par = rnorm(2), fn = lik, XXt.val = XXt.val, Uty = Uty,
#       method = "BFGS")$par)
#
# test <- c(as.numeric(-sum(XXt.val^2/v^2) + 2*sum(Uty^2*XXt.val^2/v^3)),
#           as.numeric(-sum(1/v^2) + 2*sum(Uty^2/v^3)))

fq <- function(q, kurt) {
  gamma(5/q)*gamma(1/q)/(gamma(3/q)^2) - kurt
}

fpq <- function(q) {
  (gamma(5/q)*gamma(1/q)/(q^2*gamma(3/q)^2))*(6*digamma(3/q) - digamma(1/q) - 5*digamma(5/q))
}

# Use Newton's method: https://en.wikipedia.org/wiki/Newton%27s_method
#' @export
nrq <- function(kurt, sval = 0.032, tol = 10^(-12)) { # This starting value is the lowest possible
  # Kurtosis is bounded below by 1.8, so round if needed
  kurt <- ifelse(kurt <= 1.8, 1.80001, kurt)
  # Kurtosis greater than 1.8 gives a q value of 1086.091
  # Value of fpq at q = 1086.091 is about -10^(-8), so the curve *is* pretty flat at this point
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
estRegPars <-function(y, X, delta.sq = NULL, precomp = NULL, comp.q = FALSE) {


  p <- ncol(X)
  n <- nrow(X)

  XtX <- crossprod(X)
  XtX <- crossprod(X)
  C <- cov2cor(XtX)
  V <- diag(sqrt(diag(XtX/C)))
  ceval <- eigen(C)$values
  if (!is.null(delta.sq)) {
    delta.sq <- delta.sq
  } else {
    delta.sq <- max(1 - min(ceval), 0)
  }
  D.inv <- tcrossprod(crossprod(V, (C + delta.sq*diag(p))), V)
  D <- solve(D.inv)
  # A <- diag(n) - tcrossprod(tcrossprod(X, D), X)
  DXtX <- crossprod(D, XtX)
  # XtXDDXtX <- crossprod(DXtX)
  XD <- crossprod(t(X), D)
  DXtXD <- crossprod(XD)
  # AX <- crossprod(A, X)
  # XtAAX <- crossprod(AX)
  # AA <- crossprod(A)

  b <- crossprod(D, crossprod(X, y))
  # r <- crossprod(A, y)
  #
  # e.bb <- sum(diag(XtXDDXtX))/p
  # e.be <- sum(diag(DXtXD))/p
  # e.eb <- sum(diag(XtAAX))/n
  # e.ee <- sum(diag(AA))/n
  #
  # E <- rbind(c(e.bb, e.be),
  #            c(e.eb, e.ee))
  #
  # sig.2.ests <- solve(E)%*%c(mean(b^2), mean(r^2))
  # sigma.beta.sq.hat <- sig.2.ests[1]
  # sigma.epsi.sq.hat <- sig.2.ests[2]

  vpars <- rrmmle(y = y, X = X)
  sigma.beta.sq.hat <- vpars[2]
  sigma.epsi.sq.hat <- vpars[3]

  alpha.beta <- sum(diag(crossprod(DXtX)))/p
  gamma.beta <- sum(diag(DXtX^4))/p
  omega.beta <- 3*(sum(diag(crossprod(DXtX)^2)) - sum(diag(DXtX^4)))/p

  test.stat <- (mean(b^4))/(mean(b^2)^2)

  kappa.hat <- (alpha.beta^2/gamma.beta)*(test.stat - omega.beta/alpha.beta^2)
  q.hat <- ifelse(comp.q, nrq(kappa.hat), NA)

  return(list("sigma.beta.sq.hat" = sigma.beta.sq.hat,
              "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
              "kappa.hat" = kappa.hat,
              "q.hat" = q.hat,
              "test.stat" = test.stat,
              "DXtX" = DXtX,
              "DXtXD" = DXtXD))

}


