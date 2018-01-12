HDI <- function(FUN, ...)
{
  UseMethod("HDI")
}

HDI.default <- function(FUN, lower = 0, upper = 1, level = .95, eps = 1e-3){
  
  if(!is.function(FUN)) stop("Error: 'FUN' must be a function.")
  if(length(formals(FUN)) > 1) stop("Error: 'FUN' must be a 'single-argument' function.")
  if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  x <- formals(FUN)
  fun <- function(x) FUN(x)

  posterior <- function(x) fun(x)/integrate(fun, lower, upper)[[1]]
  mode <- optimize(posterior, c(lower, upper), maximum = TRUE)[[1]]
  inverse.posterior <- function(x, side = "left") {
    target <- function(y) posterior(y) - x
    ur <- switch(side,
                 left = try(uniroot(target, interval = c(lower, mode))),
                right = try(uniroot(target, interval = c(mode, upper))))
    if(inherits(ur, "try-error")) stop("Error: You may change prior hyperparameters or limits?")
    return(ur[[1]])
  }
  areafun <- function(h) {
    i1 <- inverse.posterior(h, "left")
    i2 <- inverse.posterior(h, "right")
    return(integrate(posterior, i1, i2)[[1]])
  }
  post.area <- 1
  find.lims <- function(a) {
    ur <- uniroot(function(h) areafun(h) / post.area - a,
                  c(eps, posterior(mode) - eps))
    return(ur[[1]])
  }
  f <- find.lims(level)
  return(c(inverse.posterior(f, "left"),
           inverse.posterior(f, "right")))
}

#==================================================================================================================

hdi <- function(x, ...)
{
  UseMethod("hdi")
}

hdi.default <- function(x, y, level = .95){
  
if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
  areas <- diff(x) * .5 * (head(y, -1) + tail(y, -1))
  peak <- which.max(areas)
  range <- c(peak, peak)
  found <- areas[peak]
  while(found < level) {
    if(areas[range[1]-1] > areas[range[2]+1]) {
      range[1] <- range[1]-1
      found <- found + areas[range[1]-1]
    } else {
      range[2] <- range[2]+1
      found <- found + areas[range[2]+1]
    }
  }
  val <- x[range]
  attr(val, "indexes") <- range
  attr(val, "area") <- found
  return(val)
}

#==================================================================================================================

hdir <- function(sample, ...)
{
  UseMethod("hdir")
}

hdir.default <- function(sample, level = .95){

if(1 <= level || level <= 0) stop("Error: 'level' must be between '0' and '1'.")
if(length(sample) < 1e3) stop("Error: Insufficient sample to produce 'interval' estimates.")  
  sorted <- sort(sample)
   index <- ceiling(level*length(sorted))
       n <- length(sorted)- index
   width <- numeric(n)
  for(i in 1:n){
    width[i] <- sorted[i+ index]- sorted[i]
  }
  lower <- sorted[which.min(width)]
  upper <- sorted[which.min(width)+ index]
  return(c(lower, upper))
}

#==================================================================================================================

beta.id <- function(Low, ...)
{
  UseMethod("beta.id")
}

beta.id.default <- Vectorize(function(Low, High, Cover = NA){
  
options(warn = -1)
L <- if(is.character(Low)) as.numeric(substr(Low, 1, nchar(Low)-1)) / 100 else Low
U <- if(is.character(High)) as.numeric(substr(High, 1, nchar(High)-1)) / 100 else High
  
if(L <= 0 || U >= 1) stop("NOTE: The smallest LOWER value that you can choose is \".000001\"AND the largest UPPER value is \".999999\".")
if(L >= U) stop("Put the smaller value for Low, and the larger value for High")
  
coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100 else if(is.na(Cover)) .95 else Cover
  
p1 = (1 - coverage) / 2 
p2 = 1 - p1
  
if( p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1 ){
stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.") 
    } else {

f.beta <- function(alpha, beta, x, lower = 0, upper = 1){
p <- pbeta((x-lower)/(upper-lower), alpha, beta)
      log(p/(1-p))
    }

delta <- function(fit, actual) sum((fit-actual)^2)
 
objective <- function(theta, x, prob, ...) {
      ab <- exp(theta)
      fit <- f.beta(ab[1], ab[2], x, ...)
      return (delta(fit, prob))
    }
    
x.p <- (function(p) log(p/(1-p)))(c(p1, p2))
    
sol <- nlm(objective, log(c(1e1, 1e1)), x = c(L, U), prob = x.p, lower = 0, upper = 1, typsize = c(1, 1), 
           fscale = 1e-12, gradtol = 1e-12)
    
parm <- as.numeric(exp(sol$estimate))
    
q <- qbeta(p = c(p1, p2), parm[[1]], parm[[2]])
    
is.df <- function(a, b, sig = 3) round(a, sig) != round(b, sig)
    
if(is.df(L, q[1]) || is.df(U, q[2])){
      
stop("Error: \n\tUnable to find such a prior, make sure you have selected the correct values.")

  }else{
      
return(c(alpha = parm[[1]], beta = parm[[2]]))    
    }
  } 
})
  
#===============================================================================================

cauchy.id <- function(Low, ...)
{
  UseMethod("cauchy.id")
}
  
cauchy.id.default <- Vectorize(function(Low, High, Cover = NA){

options(warn = -1)

coverage  <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100 else if(is.na(Cover)) .95 else Cover
  
p1 = (1 - coverage) / 2
p2 = 1 - p1

if(p1 <= 0 || p2 >= 1 || Low > High || p1 > p2 || coverage >= 1) {

stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
  
} else {
  
f <- function(x) {   
    y <- c(Low, High) - qcauchy(c(p1, p2), location = x[1],  scale = x[2])
  }
  
parm <- optim(c(1, 1), function(x) sum(f(x)^2), control = list(reltol = (.Machine$double.eps)))[[1]]
}

q <- qcauchy(c(p1, p2), parm[[1]], parm[[2]])

is.df = function(a, b, sig = 4) round(a, sig) != round(b, sig)

if(is.df(Low, q[1]) || is.df(High, q[2])) {
  
stop("\n\tUnable to find such a prior, make sure you have selected the correct values")
  
} else { 
  
return(c(mode = parm[[1]], scale = parm[[2]])) 
  }
})    
  
#===============================================================================================

norm.id <- function(Low, ...)
{
  UseMethod("norm.id")
}
  
norm.id.default <- Vectorize(function(Low, High, Cover = NA){

options(warn = -1)
  
coverage <- if(is.character(Cover)) as.numeric(substr(Cover, 1, nchar(Cover)-1)) / 100 else if(is.na(Cover)) .95 else Cover
  
p1 <- (1 - coverage) / 2 
p2 <- 1 - p1
  
q <- c(Low, High)  
alpha <- c(p1, p2)
  
is.df <- function(a, b, sig = 4) (round(a, sig) != round(b, sig))
  
if( p1 <= 0 || p2 >= 1 || q[1] >= q[2] || p1 >= p2 ) {

stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
  
} else {
    
beta <- qnorm(alpha)
    
parm <- solve(cbind(1, beta), q)
    
q <- qnorm(c(p1, p2), parm[[1]], parm[[2]])
}

if(is.df(Low, q[[1]]) || is.df(High, q[[2]])) {
  
  stop("\n\tUnable to find such a prior, make sure you have selected the correct values.")
} else {
  
  return(c(mean = parm[[1]], sd = parm[[2]]))
  
  }
})
  
#===============================================================================================
 
prop.priors <- function(a, ...)
{
  UseMethod("prop.priors")
}  
  
prop.priors.default <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, scale = .1, top = 1.5, show.prior = FALSE, bottom = 1){
  
d = dist.name
pr = show.prior
is.v <- function(...) lengths(list(...)) > 1  
eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
I = eq(a, b, d, lo, hi)
a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
                           
deci <- function(x, k = 3) format(round(x, k), nsmall = k)     
                            
if(!pr){   
  Bi = round(yes)
  n = round(n)                          
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  peak = numeric(loop)
  h = list()
   
  
if(any(is.v(yes, n))) stop("Error: 'yes' & 'n' must each have a length of '1'.")  
if(yes > n) stop("Error: 'yes' cannot be larger than 'n'.")
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      h[i] = list(curve(posterior, type = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 1e3))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
      peak[i] = posterior(mode)
    }
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2, .3, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    legend("topleft", rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
    box()
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
    }
    m = scale*peak + 1:loop
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  }else{
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#==========================================================================================================================

prop.hyper <- function(a, ...)
{
  UseMethod("prop.hyper")
}

prop.hyper.default <- function(a, b, lo = 0, hi = 1, dist.name, yes = 55, n = 1e2, show.prior = FALSE, pos = 3, top = 1.01){
 
is.v <- function(...) lengths(list(...)) > 1

pr = show.prior
d = dist.name  
eq <- function(...) { lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
I = eq(a, b, d, lo, hi)
a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(a)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
                             
if(!pr){  
if(any(is.v(yes, n))) stop("Error: 'yes' & 'n' must each have a length of '1'.")  
Bi = round(yes)
n = round(n)   
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
 
 original.par = par(no.readonly = TRUE)
 on.exit(par(original.par))
    
    par(mgp = c(2.2, .3, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(0, 1), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = "Credible Interval (Proportion)", font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4", xpd = NA)
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    axis(2, at = 1:length(a), lab = deci(a), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    axis(4, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    text(par('usr')[1:2], par('usr')[4], c("A", "B"), pos = 3, cex = 1.5, xpd = NA, font = 2)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0("[", I[,1], "%", ",  ", o, "%", ",  ", I[,2], "%", "]"), cex = .75, pos = pos, xpd = NA)
   }else{
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#===================================================================================================================

ab.prop.hyper <- function(a, ...)
{
  UseMethod("ab.prop.hyper")
}

ab.prop.hyper.default <- function(a, b, lo = 0, hi = 1, dist.name, add = FALSE, 
                          yes = 55, n = 1e2, col = 1, show.prior = FALSE){
  
  is.v <- function(...) lengths(list(...)) > 1
  pr = show.prior    
  d = dist.name
  
eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
I = eq(a, b, d, lo, hi)
a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop = length(d)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
                            
if(!pr){   
if(any(is.v(yes, n))) stop("Error: 'yes' & 'n' must each have a length of '1'.")
Bi = round(yes)
n = round(n)   
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]] 
      likelihood = function(x) dbinom(Bi, n, x)
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior)
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = 0:1, xlim = c(1, loop), xlab = "Prior Parameter 'A'", xaxt = "n", yaxt = "n", ylab = "Credible Interval (Proportion)", font.lab = 2, mgp = c(2.3, .3, 0), cex.lab = 1.2)
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    axis(3, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'B'", pos = 3, cex = 1.2, xpd = NA, font = 2)
    axis(1, at = 1:length(a), lab = round(a, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){ 
    p = function(x) get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1])
    curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = "Proportion", bty = "n", font.lab = 2, lwd = 2, n = 1e3, yaxs = "i", main = bquote(Proportion*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#====================================================================================================================
     
prop.diff <- function(n, ...)
{
  UseMethod("prop.diff")
}

prop.diff.default <- function(n, yes, a = 1.2, b = a, how = c("two.one", "one.two"), level = .95, top = 1, bottom = 1, scale = .1, margin = 6){
  
    n <- round(n)
  yes <- round(yes)  
  loop <- length(n)
  
  is.s <- function(...)lengths(list(...)) < 2 
  
  if(any(yes > n)) stop("Error: 'yes' cannot be larger than 'n'.") 
  if(any(is.s(n, yes))) stop("Error: 'yes' & 'n' must each have a length of '2' or larger.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(n, yes)   
  n = I[[1]] ; yes = I[[2]] 
  
  comp <- ncol(combn(loop, 2))
  eq <- function(x) c(x, rep(rev(x)[1], comp - length(x)))
  
  if(length(a) < comp) a = eq(a)
  if(length(b) < comp) b = eq(b)
                     
  message(paste0("\n CAUTION: Check to see if you have chosen your desired ", "\"", 2*comp, "\"", " pairs of 'a' and 'b'."))
                              
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  how <- match.arg(how)
  
  delta <- switch(how,
                  one.two = function(x) x[[1]] - x[[2]], 
                  two.one = function(x) x[[2]] - x[[1]])
  
   p <- list()
for(i in 1:loop){
   p[[i]] <- rbeta(1e6, a[i] + yes[i], b[i] + (n[i] - yes[i]))
  }
  
   ps <- combn(p, 2, FUN = delta)
                  
 loop <- ncol(ps)
  
  CI <- matrix(NA, loop, 2)
 den <- list()
mode <- numeric(loop)
peak <- numeric(loop)
mean <- numeric(loop)
median <- numeric(loop)                  
  sd <- numeric(loop)
from <- numeric(loop)                  
  to <- numeric(loop)


for(i in 1:loop){
     CI[i,] <- hdir(ps[, i], level = level)
     den[i] <- list(density(ps[, i], adjust = 2, n = 1e3))
    mode[i] <- den[[i]]$x[which.max(den[[i]]$y)]
    peak[i] <- den[[i]]$y[which.max(den[[i]]$y)]
    mean[i] <- mean(ps[, i])
  median[i] <- median(ps[, i])
      sd[i] <- sd(ps[, i])
    from[i] <- mean[i] - margin *sd[i]
      to[i] <- mean[i] + margin *sd[i]
 }
  
  np <- combn(seq_along(p), 2, FUN = function(x){if(how == "one.two") paste0('p', x[1], ' - p', x[2]) else paste0('p', x[2], ' - p', x[1])})
  
  leg <- if(length(n) == 2) loop else 2
  
  plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*loop), ylab = NA, xaxt = "n", yaxt = "n", xlab = "Credible Interval (Proportion Differences)", font.lab = 2, mgp = c(2, .3, 0))
  axis(1, at = axTicks(1), labels = paste0(round(axTicks(1), 2)*1e2, "%"), mgp = c(2, .3, 0))
  abline(h = 1:loop, col = 8, lty = 3)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
  axis(2, at = 1:loop, labels = np, font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
  legend("topleft", rep(rev(paste0("beta", "(", round(a, 2), ", ", round(b, 2), ")")), leg), pch = 22, title = "Priors", pt.bg = rep(loop:1, each = leg), col = rep(loop:1, each = leg), cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
  box()
  
  for(i in 1:loop){
    polygon(x = den[[i]]$x, y = scale*den[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
  }
  
  m = scale*peak + 1:loop
  segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
  points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
  I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
  text(mode, 1:loop, paste0(I[,1], "%", "       ", o, "%", "       ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
                                                 
  rownames(CI) <- paste0(np, ":")
  colnames(CI) <- c("Lower", "Upper")
  return(data.frame(mean = mean,  mode = mode, median = median, sd = sd, CI = CI))
}     
     
#====================================================================================================================

 prop.diff.eq <- function(n1, ...)
{
  UseMethod("prop.diff.eq")
}

prop.diff.eq.default <- function(n1, n2, yes1, yes2, a1 = 1.2, b1 = 1.2, a2 = a1, b2 = b1, how = c("two.one", "one.two"), pL = -.025, pU = .025, level = .95, scale = .1){
  
  ro <- function(...){ lapply(list(...), function(x) round(x))}
  I <- ro(n1, n2, yes1, yes2)
  n1 <- I[[1]] ; n2 <- I[[2]] ; yes1 <- I[[3]] ; yes2 <- I[[4]]
  
  if(any(lengths(list(get(formalArgs(prop.diff.eq))))) > 1) stop("Error: Only 'one' comparison is allowed at a time.")
  if(yes1 > n1 || yes2 > n2) stop("Error: 'yes' cannot be larger than 'n'.")
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  p1 <- rbeta(1e6, a1 + yes1, b1 + (n1 - yes1))
  p2 <- rbeta(1e6, a2 + yes2, b2 + (n2 - yes2))
  
  how <- match.arg(how)
  
  delta <- switch(how, 
                    one.two = p1 - p2, 
                    two.one = p2 - p1) 
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(xpd = NA)
  
  d <- density(delta, adjust = 2, n = 1e5)
  
  plot(d, las = 1, type = "n", col = 0, main = NA, bty = "n", zero.line = FALSE,
       xlab = if(how == "one.two") bquote(Delta[~(p[1]-p[2])]) else bquote(Delta[~(p[2]-p[1])]), 
       cex.lab = 2, ylab = NA, axes = FALSE, yaxs = "i")
  
  axis(1, at = seq(min(d$x), max(d$x), length.out = 7), labels = paste0(deci(seq(min(d$x)*1e2, max(d$x)*1e2, length.out = 7), 2), "%"), mgp = c(2, .5, 0))
  
  polygon(x = d$x, y = scale*d$y, border = NA, col = rgb(1, 1, 0, .5)) # adjustcolor(4, .3)
  
  lines(d$x, scale*d$y, lwd = 2)     
       
  legend("topleft", c(paste0("group 1: ", "beta", "(", round(a1, 2), ", ", round(b1, 2), ")"), paste0("group 2: ", "beta", "(", round(a2, 2), ", ", round(b2, 2), ")")), title = "Priors", 
         pch = 22, col = 2, cex = .7, pt.cex = .6, pt.bg = 2, bty = "n", x.intersp = .5, title.adj = .4)
  
  mode <- d$x[which.max(d$y)]
  peak <- d$y[which.max(d$y)]*scale
  
  CI <- hdir(delta, level = level)
  segments(CI[1], 0, CI[2], 0, lend = 1, lwd = 4)
  segments(mode, 0, mode, peak, lend = 1, lty = 3)
  points(mode, 0, pch = 21, cex = 1.5, bg = "cyan")
  
  axis(side = 1, at = 0, mgp = c(3, 1.1, 0), col = 0, col.axis = "magenta", tick = FALSE, line = - 1.4, cex.axis = 1.4, font = 2)
  
  text(c(CI[1], mode, CI[2]), 0, paste0(c(deci(CI[1]*1e2, 2), deci(mode*1e2, 2), deci(CI[2]*1e2, 2)), "%"), pos = 3, 
       font = 2, col = "magenta", cex = .85)
  
  f <- approxfun(d$x, d$y, yleft = 0, yright = 0)
  
  cdf <- Vectorize(function(q){
    integrate(f, -1, q)[[1]]
  })
  
# invcdf <- function(p){
#  uniroot(function(q)cdf(q) - p, range(delta))[[1]]  # Not implemented # 
# }
  
  y1 = y2 = 1.02*peak
  x.text = (pL+pU)/2
  y.text = 1.05*peak
  low.extreme <- par('usr')[3]
  
  segments(c(pL, pU), rep(low.extreme, 2), c(pL, pU), c(y1, y2), col = 'green2', lend = 1, lty = 2)
  
  segments(pL, 0, pU, 0, col = adjustcolor(3, .5), lend = 1, lwd = 40, xpd = NA) 
       
  segments(c(pL, pU), c(y1, y2), rep(x.text, 2), rep(y.text*1.015, 2), lwd = 2, col = 'magenta')
  
  text(x.text, y.text, "Practically Equivalent to ZERO", font = 2, pos = 3, col = 'darkgreen', cex = .65, xpd = TRUE)
  
  points(c(pL, pU), c(y1, y2), pch = 21, col = 'green3', bg = 'green3', cex = 1.1)
      
  ## How much is it probable that the equivalence be true in population:

  a = cdf(pL)
  b = cdf(pU)
  
  Post.in.ROPE.Y = (b - a)
  Post.in.ROPE.X = (pU - pL) / 2
  
  BB = deci(Post.in.ROPE.Y*1e2, 2)
  
  title(main = paste0("There is ", "''", BB, "%", "''", " probability that TRUE diff. is equivalent to ZERO"), cex.main = .8)
  
  if(CI[1] > pU || CI[2] < pL) {
    
    legend("topright", "NOT Practically equivalent to \"0\" ", bty = 'n', cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    
  } else
    
    if(CI[1] > pL & CI[2] < pU) {
      
      legend("topright", "Practically equivalent to \"0\" ", bty = 'n', cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
      
    } else {
      
      legend("topright", "No decision can be made ", bty = 'n', cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    }
}             
              
#====================================================================================================================              

d.priors <- function(t, ...)
{
  UseMethod("d.priors")
}

d.priors.default <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, scale = 1, margin = 7, top = .8, show.prior = FALSE, LL = -3, UL = 3, bottom = 1, prior.left = -6, prior.right = 6){
  
  is.v <- function(...) lengths(list(...)) > 1
  d = dist.name 
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] 
  s = I[[2]] 
  d = I[[3]] 
  lo = I[[4]] 
  hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)                           
  loop = length(d) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  peak = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop) 
  h = list()
  
 if(!pr){    
  if(any(is.v(t, n1, n2))) stop("Error: 't' & 'n1' & 'n2' must each have a length of '1'.") 
  N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)   
  
  options(warn = -1)
  
  for(i in 1:loop){
    p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]        
      likelihood = function(x) dt(t, df, x*sqrt(N))
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
      sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean^2)
      from[i] = mean - margin * sd
      to[i] = mean + margin * sd
      mode[i] = optimize(posterior, c(from, to), maximum = TRUE)[[1]]
      peak[i] = posterior(mode)
      CI[i,] = HDI(posterior, LL, UL)
      h[i] = list(curve(posterior, from, to, type = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 5e2))
    }
                             
    f = peak + 1:loop
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(bottom*1, top*max(f)), ylab = NA, yaxt = "n", xlab = bquote(bold("Credible Interval "(delta))), font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    legend("topleft", rev(paste0(substring(d, 2), "(", round(m, 2), ", ", round(s, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
    box()
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop)
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
    }
    a = scale*(f-1:loop)+1:loop
    segments(mode, 1:loop, mode, a, lty = 3, xpd = NA, lend = 1)
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.1, col = 4, xpd = NA)
    I = deci(CI) ; o = deci(mode)
    text(c(CI[,1], o, CI[,2]), 1:loop, c(I[,1], o, I[,2]), pos = 3, font = 2, cex = .8, xpd = NA)
  }else{
    p = function(x) { get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
    curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }
}

#========================================================================================================================

d.hyper <- function(t, ...)
{
  UseMethod("d.hyper")
}

d.hyper.default <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, LL = -3, UL = 3, pos = 3, show.prior = FALSE, top = 1.01, margin = 6, prior.left = -6, prior.right = 6){

  is.v <- function(...) lengths(list(...)) > 1
  
  d = dist.name 
 pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]] 
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)  
  
  options(warn = -1)
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  loop = length(m)
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop)
  
if(!pr){
  if(any(is.v(t, n1, n2))) stop("Error: 't' & 'n1' & 'n2' must each have a length of '1'.")
  N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
 df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
                              
  for(i in 1:loop){
    p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]  
    likelihood = function(x) dt(t, df, x*sqrt(N))
    k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
    posterior = function(x) prior(x)*likelihood(x) / k
    mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE)[[1]]
    mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
    sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean^2)
    CI[i,] = HDI(posterior, LL, UL)
    from[i] = mean - margin * sd
    to[i] = mean + margin * sd 
    }
    
  par(mgp = c(2, .5, 0), mar = c(5.1, 4.1, 4.1, 3))   
  plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(min(from), max(to)), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xlab = bquote(bold("Credible Interval "(delta))), font.lab = 2)
  abline(h = 1:loop, col = 8, lty = 3)
  segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4")  
  points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
  axis(2, at = 1:length(m), lab = deci(m), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
  axis(4, at = 1:length(s), lab = deci(s), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
  text(par('usr')[1:2], par('usr')[4], c("M", "S"), pos = 3, cex = 1.5, xpd = NA, font = 2)
  I = deci(CI) ; o = deci(mode)
  text(mode, 1:loop, paste0("[", I[,1], ",  ", o, ",  ", I[,2], "]"), pos = pos, cex = .8, xpd = NA)
}else{
p = function(x) { get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }  
}

#===================================================================================================================

ms.d.hyper <- function(t, ...)
{
  UseMethod("ms.d.hyper")
}

ms.d.hyper.default <- function(t, n1, n2 = NA, m, s, lo = -Inf, hi = Inf, dist.name, add = FALSE, 
                      col = 1, top = 6, margin = 1.01, LL = -3, UL = 3, show.prior = FALSE, prior.left = -6, prior.right = 6){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d = dist.name 
  pr = show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I = eq(m, s, d, lo, hi)
  m = I[[1]] ; s = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
    
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)

  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  options(warn = -1)
  loop = length(m) 
  CI = matrix(NA, loop, 2)
  mode = numeric(loop)
  mean = numeric(loop)
  sd = numeric(loop)
  from = numeric(loop)
  to = numeric(loop)
                              
if(!pr){   
  if(any(is.v(t, n1, n2))) stop("Error: 't' & 'n1' & 'n2' must each have a length of '1'.")
  N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
  for(i in 1:loop){
    p = function(x) get(d[i])(x, m[i], s[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
    likelihood = function(x) dt(t, df, x*sqrt(N))
    k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
    posterior = function(x) prior(x)*likelihood(x) / k
    mode[i] = optimize(posterior, c(LL, UL), maximum = TRUE)[[1]]
    mean[i] = integrate(function(x) x*posterior(x), lo[i], hi[i])[[1]]
    sd[i] = sqrt(integrate(function(x) x^2*posterior(x), lo[i], hi[i])[[1]] - mean^2)
    CI[i,] = HDI(posterior, LL, UL)
    from[i] = mean - top * sd
    to[i] = mean + top * sd
}
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = c(min(from), max(to)), xlim = c(1, margin*loop), xlab = "Prior Parameter 'M'", xaxt = "n", ylab = bquote(bold("Credible Interval "(delta))), font.lab = 2)
    abline(v = 1:loop, col = 8, lty = 3, mgp = c(2, .5, 0))
    axis(3, at = 1:length(s), lab = round(s, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .4, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'S'", pos = 3, cex = 1, xpd = NA, font = 2)
    axis(1, at = 1:length(m), lab = round(m, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
  segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
  lines(1:loop, mode, col = col, lty = 3)
  points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){
  p = function(x){ get(d[1])(x, m[1], s[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
  curve(p, prior.left, prior.right, yaxt = "n", ylab = NA, xlab = bquote(bold("Effect Size "(delta))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(delta*" ~ "*.(if(lo[1] > -Inf || hi[1] < Inf) "truncated-")*.(substring(d[1], 2))(.(round(m[1], 2)), .(round(s[1], 2)))), mgp = c(2, .5, 0))
  }
}

#==================================================================================================================

peta.priors <- function(f, ...)
{
  UseMethod("peta.priors")
}

peta.priors.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", scale = .1, top = 1.5, show.prior = FALSE, bottom = 1){
  
  is.v <- function(...) lengths(list(...)) > 1
  
  d <- dist.name  
  pr <- show.prior
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I <- eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)                                                                                                                           

loop <- length(a)
  CI <- matrix(NA, loop, 2)
mode <- numeric(loop)
peak <- numeric(loop)
   h <- list()
                              
if(!pr){  
  
 if(any(is.v(f, N, df1, df2))) stop("Error: 'f' & 'N' & 'df1' & 'df2'  must each have a length of '1'.")  
 options(warn = -1)
                              
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]  
      likelihood = function(x) df(f, df1, df2, (x * N) / (1 - x) )
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      h[i] = list(curve(posterior, type = "n", ann = FALSE, yaxt = "n", xaxt = "n", add = i!= 1, bty = "n", n = 1e3))
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      peak[i] = posterior(mode)
      CI[i,] = HDI(posterior, 0, .9999999)
    } 
    
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = 0:1, ylim = c(bottom*1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = bquote(bold("Credible Interval"~(eta[p]^2))), font.lab = 2, mgp = c(2, .5, 0))
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, lwd = 4, col = 1:loop, xpd = NA)
    axis(2, at = 1:loop, lab = substring(d, 2), font = 2, las = 1, cex.axis = .8, tck = -.006, mgp = c(2, .3, 0))
    legend("topleft", rev(paste0(substring(d, 2), "(", round(a, 2), ", ", round(b, 2), ")")), pch = 22, title = "Priors", pt.bg = loop:1, col = loop:1, cex = .7, pt.cex = .6, bg = 0, box.col = 0, xpd = NA, x.intersp = .5, title.adj = .4)
    box()
    for(i in 1:loop){
      polygon(x = h[[i]]$x, y = scale*h[[i]]$y +i, col = adjustcolor(i, .55), border = NA, xpd = NA)
    }
    m = scale*peak + 1:loop
    segments(mode, 1:loop, mode, m, lty = 3, xpd = NA, lend = 1)  
    points(mode, 1:loop, pch = 21, bg = "cyan", cex = 1.3, col = "magenta", xpd = NA)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0(I[,1], "%", "    ", o, "%", "    ", I[,2], "%"), cex = .75, pos = 3, font = 2, xpd = NA)
  }else{
p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#===================================================================================================================

peta.hyper <- function(f, ...)
{
  UseMethod("peta.hyper")
}

peta.hyper.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", show.prior = FALSE, pos = 3, top = 1.01){
  
  is.v <- function(...) lengths(list(...)) > 1
 
  d <- dist.name
 pr <- show.prior
  
  eq <- function(...) { lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  I <- eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
loop <- length(a)
  CI <- matrix(NA, loop, 2)
mode <- numeric(loop)
                               
if(!pr){  
  if(any(is.v(f, N, df1, df2))) stop("Error: 'f' & 'N' & 'df1' & 'df2'  must each have a length of '1'.")
  options(warn = -1)
                               
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) df(f, df1, df2, (x * N) / (1 - x) )
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior, 0, .9999999)
    }

    original.par = par(no.readonly = TRUE)
    on.exit(par(original.par))
    
    par(mgp = c(2.2, .3, 0), mar = c(5.1, 4.1, 4.1, 3))   
    plot(CI[, 1:2], rep(1:loop, 2), type = "n", xlim = c(0, 1), ylim = c(1, top*loop), ylab = NA, yaxt = "n", xaxt = "n", xlab = bquote(bold("Credible Interval"~(eta[p]^2))), font.lab = 2)
    abline(h = 1:loop, col = 8, lty = 3)
    axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"))
    segments(CI[, 1], 1:loop, CI[, 2], 1:loop, lend = 1, col = "red4", xpd = NA)
    points(mode, 1:loop, pch = 21, bg = "red4", cex = .8, col = "red4", xpd = NA)
    axis(2, at = 1:length(a), lab = deci(a), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    axis(4, at = 1:length(b), lab = deci(b), font = 2, las = 1, cex.axis = .8, tick = FALSE, mgp = c(2, .2, 0))
    text(par('usr')[1:2], par('usr')[4], c("A", "B"), pos = 3, cex = 1.5, xpd = NA, font = 2)
    I = deci(CI*1e2 , 2); o = deci(mode*1e2, 2)
    text(mode, 1:loop, paste0("[", I[,1], "%", ",  ", o, "%", ",  ", I[,2], "%", "]"), cex = .75, pos = pos, xpd = NA)
  }else{
p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}

#===================================================================================================================

ab.peta.hyper <- function(f, ...)
{
  UseMethod("ab.peta.hyper")
}

ab.peta.hyper.default <- function(f, N, df1, df2, a = 1.2, b = 1.2, lo = 0, hi = 1, dist.name = "dbeta", add = FALSE, 
                          col = 1, show.prior = FALSE){
   
  is.v <- function(...) lengths(list(...)) > 1
   
  d <- dist.name
  pr <- show.prior    
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
   I <- eq(a, b, d, lo, hi)
  a = I[[1]] ; b = I[[2]] ; d = I[[3]] ; lo = I[[4]] ; hi = I[[5]]
  
  if(is.v(a) & pr || is.v(b) & pr) message("\tNote: You can see only '1 prior' at a time.")
  if(add & pr) message("\tNote: 'add' only works for overlying 'Credible Intervals' to compare them.")
  
  deci <- function(x, k = 3) format(round(x, k), nsmall = k)
  
  loop <- length(a)
    CI <- matrix(NA, loop, 2)
  mode <- numeric(loop)
      
 options(warn = -1)
         
 if(!pr){    
  if(any(is.v(f, N, df1, df2))) stop("Error: 'f' & 'N' & 'df1' & 'df2'  must each have a length of '1'.")                       
  for(i in 1:loop){
    p = function(x) get(d[i])(x, a[i], b[i])*as.integer(x >= lo[i])*as.integer(x <= hi[i])
    prior = function(x) p(x)/integrate(p, lo[i], hi[i])[[1]]
      likelihood = function(x) df(f, df1, df2, (x * N) / (1 - x) )
      k = integrate(function(x) prior(x)*likelihood(x), lo[i], hi[i])[[1]]
      posterior = function(x) prior(x)*likelihood(x) / k
      mode[i] = optimize(posterior, c(lo[i], hi[i]), maximum = TRUE)[[1]]
      CI[i,] = HDI(posterior, 0, .9999999)
    }
  }
  
  if(!add & !pr){
    plot(rep(1:loop, 2), CI[, 1:2], type = "n", ylim = 0:1, xlim = c(1, loop), xlab = "Prior Parameter 'A'", xaxt = "n", yaxt = "n", ylab = bquote(bold("Credible Interval"~(eta[p]^2))), font.lab = 2, mgp = c(2.3, .3, 0), cex.lab = 1.2)
    abline(v = 1:loop, col = 8, lty = 3)
    axis(2, at = axTicks(2), lab = paste0(axTicks(2)*1e2, "%"), mgp = c(2, .4, 0), las = 1)
    axis(3, at = 1:length(b), lab = round(b, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .2, 0))
    text(mean(par('usr')[1:2]), 1.06*par('usr')[4], "Prior Parameter 'B'", pos = 3, cex = 1.2, xpd = NA, font = 2)
    axis(1, at = 1:length(a), lab = round(a, 3), font = 2, las = 1, cex.axis = .8, mgp = c(2, .3, 0))
  }
  
  if(!pr){
    segments(1:loop, CI[, 1], 1:loop, CI[, 2], lend = 1, col = col)  
    lines(1:loop, mode, col = col, lty = 3)
    points(1:loop, mode, pch = 21, bg = col, cex = .8, col = col, xpd = NA)
  }
  
  if(!add & pr){
p = function(x) { get(d[1])(x, a[1], b[1])*as.integer(x >= lo[1])*as.integer(x <= hi[1]) }
curve(p, 0, 1, yaxt = "n", xaxt = "n", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), bty = "n", font.lab = 2, lwd = 2, n = 1e3, main = bquote(eta[p]^2*" ~ "*.(if(lo[1] > 0 || hi[1] < 1) "truncated-")*.(substring(d[1], 2))(.(round(a[1], 2)), .(round(b[1], 2)))))
axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .4, 0))
  }
}


#=================================================================================================================

prop.update <- function(n, ...)
{
  UseMethod("prop.update")
}

prop.update.default <- function(n = 100, yes = 55, top = 5, scale = .1, lo = 0, hi = 1, a = 1.5, b = 1.5, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5){

pri <- show.prior
s <- yes  
d <- dist.name 
is.v <- function(...) lengths(list(...)) > 1
if(any(is.v(a, b, d))) stop("Error: Choose only 'one' prior knowledge base at a time.")  
if(d == "dunif" & a == 0 & b == 1) {d = 'dbeta'; a <- b <- 1.0000001}
if(d == "dbeta" & a == 1 & b == 1) a <- b <- 1.0000001;    
if(tol < 1e4) stop("'tol' must be '10,000' or larger.")

eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
I <- eq(n, s) ; n <- I[[1]] ; s <- I[[2]]
loop <- length(n) 

  props <- seq(0, 1, 1/tol)
  prx <- get(d)(props, a, b)*as.integer(props >= lo)*as.integer(props <= hi)
  pr <- tol * prx / sum(prx)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
   
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~props, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = NA, xlab = "Proportion", font.lab = 2, axes = FALSE, mgp = c(2, .4, 0),main = if(pri) bquote(Proportion*" ~ "*.(if(lo > 0 || hi < 1) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))

if(!pri){  
  abline(h = 1:loop+1, col = 8, lty = 3)
  axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
}  
  polygon(x = c(lo, props, hi), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  I = hdi(x = props, y = pr, level = level)
  
  m = props[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(m, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(m, 0, m, y, lty = 3)
  q <- deci(I*1e2, 2) 
  o <- deci(m*1e2, 2)
  text(c(I[1], m, I[2]), 0, c(paste0(q[1], "%"), paste0(o, "%"), paste0(q[2], "%") ), pos = 3, cex = .8, font = 2, xpd = NA)

if(!pri){
  for(i in 1:loop) {
    ps <- dbinom(s[i], n[i], props) * pr
    ps <- tol * ps / sum(ps)
polygon(y = scale*ps+i+1, x = props, col = adjustcolor(i+1, .5), border = NA, xpd = NA)
I = hdi(x = props, y = ps, level = level)
m = props[which.max(ps)]
q = deci(I*1e2 , 2); o = deci(m*1e2, 2)
y = ps[which.max(ps)]*scale + (i+1)
segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
segments(m, i+1, m, y, lty = 3, xpd = NA)
text(m, i+1, paste0(q[1], "%", "     ", o, "%", "     ", q[2], "%"), pos = 3, cex = .7, font = 2, xpd = NA)
points(m, i+1, pch = 21, bg = "cyan", col = "magenta")

    pr <- ps
    }
  }
}

#=======================================================================================================================

d.update <- function(t, ...)
{
  UseMethod("d.update")
}

d.update.default <- function(t, n1, n2 = NA, top = 5, scale = .1, m = 0, s = 1, dist.name, prior.scale = 1, level = .95, show.prior = FALSE, lo = -2, hi = 2, tol = 1e4, margin = hi){
  
  pri <- show.prior
  d <- dist.name
  if(is.infinite(lo)) lo <- -6
  if(is.infinite(hi)) hi <- 6
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(m, s, d))) stop("Choose only 'one' prior knowledge base at a time.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(t, n1, n2) 
  t <- I[[1]]  
  n1 <- I[[2]]  
  n2 <- I[[3]] 
  loop <- length(t) 
  
  ds <- seq(-6, 6, 1/tol)
  prx <- get(d)(ds, m, s)*as.integer(ds >= lo)*as.integer(ds <= hi)
  pr <- tol * prx / sum(prx)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~ds, ylim = c(0, top*loop), xlim = c(-margin, margin), type = "n", xaxs = "i", yaxs = "i", ylab = NA, xlab = bquote(bold("Effect Size"~ (delta))), font.lab = 2, mgp = c(2, .4, 0),main = if(pri) bquote("Effect Size"*" ~ "*.(if(lo > -Inf || hi < Inf) "truncated-")*.(substring(d, 2))(.(round(m, 2)), .(round(s, 2)))) else NA, yaxt = "n", bty = "n")
  
  if(!pri){  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
  }  

  polygon(x = c(-margin, ds, margin), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  # I = hdi(x = ds, y = pr, level = level)
  
  if(d != "dunif"){
  mode = ds[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  # segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(mode, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(mode, 0, mode, y, lty = 3)
  # text(c(.85*I[1], mode, I[2]), 0, paste0(round(c(I[1], mode, I[2]), 3)), pos = 3, cex = .8, font = 2, xpd = NA)
  }
  
  N = ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
  df = ifelse(is.na(n2), n1 - 1, n1 + n2 - 2) 
  
  options(warn = -1)
  if(!pri){
    for(i in 1:loop) {
     
      ps <- dt(t[i], df[i], ds*sqrt(N[i])) * pr
      ps <- tol * ps / sum(ps)
      polygon(y = scale*ps+i+1, x = ds, col = adjustcolor(i+1, .5), border = NA, xpd = NA)
      I = hdi(x = ds, y = ps, level = level)
      mode = ds[which.max(ps)]
      q = deci(I, 3); o = deci(mode, 3)
      y = ps[which.max(ps)]*scale + (i+1)
      segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
      segments(mode, i+1, mode, y, lty = 3, xpd = NA)
      text(mode, i+1, paste0(q[1], "     ", o, "     ", q[2]), pos = 3, cex = .7, font = 2, xpd = NA)
      points(mode, i+1, pch = 21, bg = "cyan", col = "magenta")
      
      pr <- ps
    }
  }
}

#==================================================================================================================

peta.update <- function(f, ...)
{
  UseMethod("peta.update")
}

peta.update.default <- function(f, N, df1, df2, top = 5, scale = .1, a = 2, b = 2, lo = 0, hi = 1, dist.name = "dbeta", prior.scale = 1, level = .95, show.prior = FALSE, tol = 1e5){
  
  pri <- show.prior
  d <- dist.name
  if(hi == 1) hi <- .9999999 ;
  if(tol < 1e4) stop("'tol' must be '10,000' or larger.")
  is.v <- function(...) lengths(list(...)) > 1
  if(any(is.v(a, b, d))) stop("Error: Choose only 'one' prior knowledge base at a time.")
  
  eq <- function(...){ lapply(list(...), function(x) c(x, rep(rev(x)[1], max(lengths(list(...))) - length(x)))) }
  deci <- function(x, k = 3) format(round(x, k), nsmall = k) 
  I <- eq(f, N, df1, df2) ; f <- I[[1]] ; N <- I[[2]] ; df1 <- I[[3]] ; df2 <- I[[4]] ;
  loop <- length(f) 
  
  peta <- seq(0, .9999999, 1/tol)
  prx <- get(d)(peta, a, b)*as.integer(peta >= lo)*as.integer(peta <= hi)
  pr <- tol * prx / sum(prx)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(5, 6.8, 4, 2))
  plot(pr~peta, ylim = c(0, top*loop), type = "n", yaxs = "i", ylab = NA, xlab = bquote(bold("Partial Eta.Sq"~(eta[p]^2))), font.lab = 2, axes = FALSE, mgp = c(2, .4, 0),main = if(pri) bquote(eta[p]^2*" ~ "*.(if(lo > 0 || hi < .9999999) "truncated-")*.(substring(d, 2))(.(round(a, 2)), .(round(b, 2)))) else NA)
  axis(1, at = axTicks(1), lab = paste0(axTicks(1)*1e2, "%"), mgp = c(2, .3, 0))
  
  if(!pri){  
    abline(h = 1:loop+1, col = 8, lty = 3)
    axis(2, at = 0:loop+1, lab = c("Base knowledge", paste0("Study ", 1:loop)), las = 1, font = 2, cex.axis = .9, mgp = c(2, .2, 0), tick = FALSE, xpd = NA)
  }  
  polygon(x = c(lo, peta, hi), y = prior.scale*c(0, pr, 0), col = adjustcolor(8, .8))
  
  I = hdi(x = peta, y = pr, level = level)
  
  m = peta[which.max(pr)]
  y = prior.scale*(pr[which.max(pr)])
  segments(I[1], 0, I[2], 0, lend = 1, lwd = 4, xpd = NA)
  points(m, 0, pch = 19, xpd = NA, cex = 1.4)
  segments(m, 0, m, y, lty = 3)
  text(c(.85*I[1], m, I[2]), 0, paste0(round(c(I[1], m, I[2])*1e2, 4), "%"), pos = 3, cex = .8, font = 2, xpd = NA)
  
  if(!pri){
    for(i in 1:loop) {
      
      ps <- df(f[i], df1[i], df2[i], (peta * N[i]) / (1 - peta) ) * pr
      ps <- tol * ps / sum(ps)
      m = peta[which.max(ps)]
      polygon(y = c(i+1, scale*ps+i+1, i+1), x = c(lo, peta, hi), col = adjustcolor(i+1, .5), border = NA, xpd = NA)
      I = hdi(x = peta, y = ps, level = level)
      
      q = deci(I*1e2 , 2); 
      o = deci(m*1e2, 2)
      y = ps[which.max(ps)]*scale + (i+1)
      segments(I[1], i+1, I[2], i+1, lend = 1, lwd = 3, col = i +1)
      segments(m, i+1, m, y, lty = 3, xpd = NA)
      text(m, i+1, paste0(q[1], "%", "     ", o, "%", "     ", q[2], "%"), pos = 3, cex = .7, font = 2)
      points(m, i+1, pch = 21, bg = "cyan", col = "magenta")
      
      pr <- ps
    }
  }
}

#===================================================================================================================

d.eq.test <- function(t, ...)
{
  UseMethod("d.eq.test")
}

d.eq.test.default <- function(t, n1, n2 = NA, m, s, dist.name, dL = -.1, dU = .1, lo = -Inf, hi = Inf){
  
  d <- dist.name
  
  if(any(lengths(list(get(formalArgs(d.eq.test))))) > 1) stop("Error: Only 'one' equivalence testing at a time is allowed.")
  if(dL >= dU) stop("Your Upper value must be larger than your Lower value")
  
  if(abs(dL) != abs(dU)) message("\n\tYou have an \"Unequal Equivalence Bound\", thus we can't provide an extra\n\t function showing the effect of choosing various Unequal bounds.")
  
  decimal <- function(x, k){   
    if(is.character(x)){
      return(x)
    }else{
    format(round(x, k), nsmall = k, scientific =
    ifelse(x >= 1e5 || x <= -1e5 || x <= 1e-5 & x >= -1e-5, TRUE, FALSE) )
       }
    }
 
  options(warn = -1)
  
   N <- ifelse(is.na(n2), n1, (n1 * n2) / (n1 + n2))
  df <- ifelse(is.na(n2), n1 - 1, n1 + n2 - 2)
  
         p <- function(x) get(d)(x, m, s)*as.integer(x >= lo)*as.integer(x <= hi)
     prior <- function(x) p(x)/integrate(p, lo, hi)[[1]]  
likelihood <- function(x) dt(t, df, x*sqrt(N))
         k <- integrate(function(x) prior(x)*likelihood(x), lo, hi)[[1]]
 posterior <- function(x) prior(x)*likelihood(x)/k
  
  mean <- integrate(function(x) x*posterior(x), lo, hi)[[1]]
    sd <- sqrt(integrate(function(x) x^2*posterior(x), lo, hi)[[1]] - mean^2)
    
  x.min.1 <- mean - 9 * sd
  x.max.1 <- mean + 9 * sd
  
  ## The dL and dU may be different from x.min.1 and x.max.1 respectively, if so, adjust accordingly.
  x.min <- if(dL < x.min.1) { 1.05*dL } else { x.min.1 }
  x.max <- if(dU > x.max.1) { 1.05*dU } else { x.max.1 }
  
  CI <- HDI(posterior, x.min, x.max)
  
  cdf <- Vectorize(function(q){
    integrate(posterior, lo, q)[[1]]
  })
  
#  inv.cdf <- Vectorize(function(p){
#    uniroot(function(q)cdf(q) - p, c(x.min, x.max))[[1]]  # Not implemented #
#  })
  
  mode <- optimize(posterior, c(x.min, x.max), maximum = TRUE)[[1]]
  peak <- posterior(mode)
  
  original.par = par(no.readonly = TRUE)
  on.exit(par(original.par))
  
  par(mar = c(.1, 4.1, 3.1, 2.1), mfcol = c(2, 1))
  
  h = curve(posterior, from = x.min, to = x.max, las = 1, type = "n",
             xlab = NA, ylab = NA, bty = "n", ylim = c(0, 1.1*peak), 
             xaxt = "n", yaxt = "n", mgp = c(2, .5, 0), n = 1e3)
  
  X <- h$x >= CI[1] &  h$x <= CI[2]
  
  low.extreme <- par('usr')[3]
  
  polygon(c(CI[1], h$x[X], CI[2]), c(low.extreme, h$y[X], low.extreme), col = rgb(1, 1, 0, .5), border = NA)
  
  segments(mode, low.extreme, mode, peak, lty = 3)
  
  text(mode, peak/2, decimal(mode, 2), srt = 90, pos = 3, font = 2)
  
  lines(h, lwd = 2)
  
  segments(CI[1], low.extreme, CI[2], low.extreme, col = 2, lend = 1, lwd = 40)
  
  segments(dL, low.extreme, dU, low.extreme, col = adjustcolor(3, .5), lend = 1, lwd = 40)
  
  points(mode, low.extreme/5, pch = 21, col = 0, bg = 0, cex = 1.5)
  
  axis(side = 1, at = decimal(seq(x.min, x.max, length.out = 7), 2), mgp = c(2, .5, 0))
  axis(side = 1, at = 0, mgp = c(3, 1.1, 0), col = 0, col.axis = "magenta", tick = FALSE, line = - 1.4, cex.axis = 1.4, font = 2)
  
  mtext(side = 1, bquote(bold("Population Effect Size"~(delta))), line = 2, cex = .95)
  
  y1 = y2 = 1.02*peak
  x.text = (dL+dU)/2
  y.text = 1.05*peak
  
  segments(c(dL, dU), rep(low.extreme, 2), c(dL, dU), c(y1, y2), col = 'green2', lend = 1, lty = 2)
  
  segments(c(dL, dU), c(y1, y2), rep(x.text, 2), rep(y.text*1.023, 2), lwd = 2, col = 'magenta')
  
  text(x.text, y.text, "Practically Equivalent to ZERO", font = 2, pos = 3, col = 'darkgreen', cex = .65, xpd = TRUE)
  
  points(c(dL, dU), c(y1, y2), pch = 21, col = 'green3', bg = 'green3', cex = 1.1)
  
  ## How much is it probable that the equivalence be true in population:
  
  a = cdf(dL)
  b = cdf(dU)
  
  Post.in.ROPE.Y = (b - a)
  Post.in.ROPE.X = (dU - dL) / 2
  
  BB = decimal(Post.in.ROPE.Y*1e2, 2)
  
  title(main = paste0("There is ", "''", BB, "%", "''", " probability that TRUE effect size is equivalent to ZERO"), cex.main = .8)
  
  legend("topleft", legend = paste0("95% HDI: [",decimal(CI[1], 2), ", ", decimal(CI[2], 2),"]"),
         bty = "n", inset = c(-.035,.1), text.font = 2, text.col = 'red4', cex = .8)
  
  if(CI[1] > dU || CI[2] < dL) {
    
    legend("topright", "NOT Practically equivalent to \"0\" ", bty = 'n', inset = c(-.01, .1), cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    
  } else
    
    if(CI[1] > dL & CI[2] < dU) {
      
      legend("topright", "Practically equivalent to \"0\" ", bty = 'n', inset = c(-.01, .1), cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
      
    } else  {
      
      legend("topright", "No decision can be made ", bty = 'n', inset = c(-.01, .1), cex = .75, text.font = 4, text.col = 'magenta2', title = "Decision:")
    }
  
  #########################################################################
  ## How choice of ROPE can affect porortion of posterior that ROPE covers:
  #########################################################################
  
  par(mar = c(3.1, 4.1, 6.1, 2.1), mgp = c(2.5, .5, 0))
  
  eq.low = ifelse(abs(dL) <= .3, 4, 2)*( - ((dU - dL) / 2) )
  eq.upp = ifelse(abs(dL) <= .3, 4, 2)*(   ((dU - dL) / 2) )

  L = seq(eq.low, 0, length.out = 1e2)
  U = seq(eq.upp, 0, length.out = 1e2)
    
  aa = cdf(L)
  bb = cdf(U)
  
  Eq = (bb - aa)  # porortion of posterior that ROPE covers
  half = (U - L)/2
  
  plot(half, Eq, type = ifelse(abs(dL) == abs(dU), 'l' ,'n'), lwd = 3, col = 'red4', axes = FALSE,
       xlab = NA, ylab = paste0("%",'Posterior in ROPE'), font.lab = 2, cex.lab = .8)
  
  mtext(side = 1, "Half of ROPE", font = 2, line = 1.5, cex = .9)
  
  axis(1, at = decimal(seq(0, eq.upp[1], length.out = 7), 2), las = 1)
  axis(2, at = seq(0, Eq[1], length.out = 5),
       labels = paste0(1e2*round(seq(0, Eq[1], length.out = 5),
                                2), "%"), las = 1)
  
  pars = par('usr')
  
  rect(pars[1], pars[3], pars[2], pars[4], col = adjustcolor("grey", .1), border = NA)
    
  rect(pars[1], pars[3], Post.in.ROPE.X, Post.in.ROPE.Y,
       col = adjustcolor("yellow",
                         alpha.f = ifelse(Post.in.ROPE.Y  <= .2, .1,
                                          ifelse(Post.in.ROPE.Y > .2 & Post.in.ROPE.Y <= .3, .15,
                                                 ifelse(Post.in.ROPE.Y > .3 & Post.in.ROPE.Y <= .4, .2, .3)))),
       lty = 2 )
    
  points(Post.in.ROPE.X, Post.in.ROPE.Y,
         pch = 21, cex = 2, bg = 'green')
  
  box()  
}
