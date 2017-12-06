HDI <- function(Posterior, domain = c(0, 1), level = .95, eps = 1e-3) {
  
  lower = min(domain) ; upper = max(domain)
  posterior = function(x) Posterior(x)/integrate(Posterior, lower, upper)[[1]]
  mode = optimize(posterior, interval = domain, maximum = TRUE, tol = 1e-20)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("inverse.posterior failed: You may change prior specification or extend limit?")
    return(ur$root)
  }
  areafun <- function(h) {
    i1 <- inverse.posterior(h, "left")
    i2 <- inverse.posterior(h, "right")
    return(integrate(posterior, i1, i2)[[1]])
  }
  post.area <- integrate(posterior, lower, upper)[[1]]
  if(post.area<level) stop("limits don't encompass desired area: a =", round(post.area, 3))
  find.lims <- function(a) {
    ur <- uniroot(function(h) areafun(h) / post.area - a,
                  c(eps, posterior(mode) - eps))
    return(ur$root)
  }
  f <- find.lims(level)
  return(c(inverse.posterior(f, "left"),
           inverse.posterior(f, "right")))
}
