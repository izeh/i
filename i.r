HDI <- function(FUN, lower = 0, upper = 1, level = .95, eps = 1e-3){
  
  if(!is.function(FUN)) stop("Error: 'FUN' must be a function.")
  if(length(formals(FUN)) > 1) stop("Error: 'FUN' must be a 'single-argument' function.")
  x <- formals(FUN)
  fun <- function(x) FUN(x)

  lower = min(lower, upper) ; upper = max(lower, upper)
  posterior = function(x) fun(x)/integrate(fun, lower, upper)[[1]]
  mode = optimize(posterior, interval = c(lower, upper), maximum = TRUE, tol = 1e-20)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("inverse.posterior failed: You may change prior hyperparameters or extend limit?")
    return(ur$root)
  }
  areafun <- function(h) {
    i1 <- inverse.posterior(h, "left")
    i2 <- inverse.posterior(h, "right")
    return(integrate(posterior, i1, i2)[[1]])
  }
  post.area <- 1
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
