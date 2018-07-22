
need <- "gsl"
have <- need %in% rownames(installed.packages())
if(!have){ install.packages( need[!have] ) }

options(warn = -1)
suppressMessages({ 
  library("gsl")
})


peta2F <- function(peta, df1, df2) (peta / df1) / ((1 - peta)/df2)

peta2ncp <- function(peta, N)  (peta*N) / (1 - peta) 

ncp2peta <- function(ncp, N) ncp / (ncp + N)


peta.ci <- function(peta, f = NA, df1, df2, N, conf.level = .9, digits = 20)
{
  UseMethod("peta.ci")
} 

peta.ci.default <- function(peta, f = NA, df1, df2, N, conf.level = .9, digits = 9){
  
  ci <- Vectorize(function(peta, f, N, df1, df2, conf.level){
    
    q <- ifelse(is.na(f), peta2F(peta, df1, df2), f) 
    
    alpha <- (1 - conf.level)/2
    
    u <- function (ncp, alpha, q, df1, df2) {
      suppressWarnings(pf(q = q, df1 = df1, df2 = df2, ncp, lower.tail = FALSE)) - alpha
    }
    
    g <- try(uniroot(u, c(0, q+1e7), alpha = alpha, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], silent = TRUE)
    if(inherits(g, "try-error")) g <- 0
    h <- try(uniroot(u, c(0, q+1e7), alpha = 1-alpha, q = q, df1 = df1, df2 = df2, extendInt = "yes")[[1]], silent = TRUE)
    if(inherits(h, "try-error")) h <- 0
    I <- c(g, h)
    
    I <- I / (I + N)
    
    P.eta.sq <- if(is.na(f)) peta else F2peta(f, df1, df2)
    
    return(c(P.eta.sq = P.eta.sq, lower = I[1], upper = I[2], conf.level = conf.level, ncp = peta2ncp(P.eta.sq, N), F.value = q))
  })
  
  peta <- if(missing(peta)) NA else peta
  
  round(data.frame(t(ci(peta = peta, f = f, N = N, df1 = df1, df2 = df2, conf.level = conf.level))), digits = digits)
}


exp.pov <- function(P2, K, N, regress = TRUE)
{
  K <- if(regress) K else K + 1 
  expect <- 1 - ((N - K - 1)/(N - 1)) * (1 - P2) * gsl::hyperg_2F1(1, 1, (N + 1)/2, P2)
  max(0, expect)
}


"%inn%" <- function(x = 3.5, interval = c(3, 5)){
  
  r <- range(interval, na.rm = TRUE)
  
  x >= r[1] & x <= r[2] 
}


root <- function(pov = .6, df1 = 3, df2 = 108, N = 100, conf.level = .95, show = FALSE, ...){
  
  f <- function(x){ 
    
    ci <- peta.ci(peta = x, df1 = df1, df2 = df2, N = N, conf.level = conf.level, digits = 1e2)
    
    abs(ci$upper - ci$lower)
  }
  
  m <- optimize(f, 0:1, maximum = TRUE)[[1]]
  
  est <- uniroot(function(x) f(pov) - f(x), if(pov >= m) c(0, m) else c(m, 1))[[1]]
  
  if(show) curve(f, panel.f = abline(v = c(pov, est), h = f(pov), col = 2, lty = c(2, 1, 1)), ...) 
  
  list(m = m, est = est)
}


plan.f.ci <- function(H2 = .2, design = 2 * 2, n.level = 2, n.covar = 0, conf.level = .9, width = .2, regress = FALSE,  pair.design = 0, assure = .99){
  
  UseMethod("plan.f.ci")
}
                 
                 
plan.f.ci.default <- function(H2 = .2, design = 2 * 2, n.level = 2, n.covar = 0, conf.level = .9, width = .2, regress = FALSE,  pair.design = 0, assure = .99){
  
  if(any(conf.level >= 1) || any(conf.level <= 0) || any(assure >= 1) || any(assure <= 0)) stop("'conf.level' and 'assure' must be between '0' and '1'.", call. = FALSE)
  peta <- H2
  G <- Vectorize(function(peta, conf.level, width, assure, design, n.level, n.covar, regress, pair.design){
    
    n.f <- function(peta, conf.level, width, assure, design, n.level, n.covar, regress, pair.design){
      
      alpha <- (1 - conf.level)/2
      if(regress){ n.level <- n.level + 1 ; design <- n.level }
      if(pair.design != 0) design <- 2
      df1 <- n.level - 1
      if(n.covar < 0) n.covar <- 0
      options(warn = -1)
      
      f <- function(alpha, q, df1, df2, ncp){
        alpha - suppressWarnings(pf(peta2F(peta, df1, df2), df1, df2, ncp, lower.tail = FALSE))
      }
      
      pbase <- function(df2){      
        
        b <- sapply(c(alpha, 1 - alpha), function(x) 
          tryCatch(uniroot(f, c(0, 1e7), alpha = x, q = q, df1 = df1, df2 = df2)[[1]], error = function(e) NA))
        if(any(is.na(b))) b <- c(1, 1e4)     
        ncp2peta(b, df2 + design)
      }
      
      m <- function(df2, width){
        abs(diff(pbase(df2))) - width
      }
      
      df2 <- uniroot(m, c(0, 1e3), width = width, extendInt = "yes")[[1]]
      
      df2 <- if(regress) ceiling(df2) else ceiling(df2 - n.covar)
      
      N <- ceiling(df2 + design)
      bal <- ceiling(N/design) * design
      if(pair.design != 0){ N <- pair.design * (bal/2) ; message("\nNote: You are doing reseach planning for accurate 'pairwise' comparisons.") }
      N <- if(!regress & design != 0 & N %% design != 0) bal else N
      n.covar <- if(n.covar == 0) NA else n.covar
      n.level <- if(regress) n.level-1 else n.level
      design <- if(regress) n.level else design
      df1 <- if(regress) n.level else df1
      
      list(peta = peta, total.N = N, width = width, n.level = n.level, conf.level = conf.level, assure = assure, df1 = df1, df2 = df2, design = design)
    }
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar, pair.design = pair.design)  
    
    peta <- exp.pov(P2 = n$peta, K = n$design, N = n$total.N, regress = regress)
    
    n <- n.f(peta = peta, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar, pair.design = pair.design)
    
    peta.max <- root(pov = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = conf.level)$m
    
    a <- peta.ci(peta = peta, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = assure - (1 - assure))
    
    nLU <- sapply(c(a$lower, a$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar, pair.design = pair.design)$total.N)
    
    NN1 <- max(nLU, na.rm = TRUE)
    
    b <- peta.ci(peta = peta.max, df1 = n$df1, df2 = n$df2, N = n$total.N, conf.level = 1 - assure)
    
    nLU <- sapply(c(b$lower, b$upper), function(x) n.f(peta = x, width = width, assure = assure, n.level = n.level, regress = regress, conf.level = conf.level, design = design, n.covar = n.covar, pair.design = pair.design)$total.N)
    
    NN2 <- max(nLU, na.rm = TRUE)
    
    NN3 <- if(!(peta.max %inn% c(a$lower, a$upper))) NN1 else max(NN1, NN2)
    
    return(c(peta = peta, total.N = NN3, width = width, n.level = n.level, design = design, conf.level = conf.level, N1 = NN1, N2 = NN2))
    
  })
  
  data.frame(t(G(peta = peta, conf.level = conf.level, width = width, design = design, n.level = n.level, n.covar = n.covar, pair.design = pair.design, assure = assure, regress = regress)), regress = regress, row.names = NULL)
}



