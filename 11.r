
pre.post.cont <- function(n1 = 12, n2 = 12, min.score = 0, max.score = 25, subjects = TRUE, 
                  descriptives = TRUE, correlation = .7, effect.size = 1, digits = 6, ...){
  
  decimal <- function(x, k){
    if(is.character(x)){ x 
    }else{
      format(round(x, k), nsmall = k, scientific =
               ifelse(x >= 1e5 || x <= -1e5 || x <= 1e-5 & x >= -1e-5, TRUE, FALSE) )
    }
  }  
  
  if(min.score >= max.score){
    stop("\n\tYour \"min.score\" must be smaller than your \"max.score\".")  }
  
  beta = qnorm(c(1e-10, .9999999999))
  q = c(min.score, max.score)
  
  mu.sigma = solve(cbind(1L, beta), q)
  
  mean = mu.sigma[[1]]
  sd = mu.sigma[[2]]
  
  coeff = effect.size*sd
  
  aa = mean + .5*mean
  bb = mean + .3*mean
  cc = aa - coeff
  
  mean.g2 = min(bb, cc)
  mean.g1 = mean.g2 + coeff
  
  TRUE.d = (mean.g1 - mean.g2) / sd      
  
  cor.pop = correlation
  
  mu <- c(0, 0)
  cov.pop <- matrix(c(1, cor.pop, cor.pop, 1), nrow = 2)
  
  need <- "MASS"
  have <- need %in% rownames(installed.packages())
  if(any(!have)){ install.packages( need[!have] ) }
  
  suppressMessages({ 
    library("MASS")
  })
  mvnorm.mat <- mvrnorm(n1, Sigma = cov.pop, mu = mu)
  
  a <- mvnorm.mat[ , 1] * sd + mean.g1
  b <- mvnorm.mat[ , 2] * sd + mean.g2
  
  y1 = c(a - b)
  
  mvnorm.mat <- mvrnorm(n2, Sigma = cov.pop, mu = mu)
  
  a <- mvnorm.mat[ , 1] * sd + mean.g2
  b <- mvnorm.mat[ , 2] * sd + mean.g2
  
  y2 = c(a - b)
  
  y = c(y1, y2)
  
  groups = factor(rep(1:2, c(n1, n2)), labels = c("Treatment", "Control"))
  
  mean.g1 = mean(y[groups == "Treatment"])  
  mean.g2 = mean(y[groups == "Control"])    
  
  sd.g1 = sd(y[groups == "Treatment"])
  sd.g2 = sd(y[groups == "Control"])
  
  groups.for.t = factor(rep(1:2, c(n1, n2)))
  
  test = t.test(y ~ groups.for.t, var.equal = TRUE)
  
  t.value = unname(test$statistic) ; p.value = test$p.value
  
  Cohend = t.value / sqrt((n1*n2)/(n1+n2)) 
  
  lab1 = if(n1 < 10 || n2 < 10) paste0("subj #", rev(1L:n1)) else c(paste0("subj #", rev(n1)[1]), paste0(rep(".", n1 - 2)), paste0("subj #", 1L))
  lab2 = if(n1 < 10 || n2 < 10) paste0("subj #", rev(1L:n2)) else c(paste0("subj #", rev(n2)[1]), paste0(rep(".", n2 - 2)), paste0("subj #", 1L))
  
  if(subjects) {
    graphics.off()
    par(font.lab = 2, mar = c(4.2, 3, 2, 4), xaxt = "n", ...)
    dotchart(y, groups = groups, color = c(4, 2)[groups], 
             font = 2, pch = 19, gcolor = c(4, 2), xlab = "Participants' Gain Scores",
             pt.cex = ifelse(n1 <= 20 || n2 <= 20, 1.5, .8), labels = c(lab1, lab2), main = NA,
             cex.main = 2) 
  } else {
    graphics.off()
    par(font.lab = 2, mar = c(4.2, 3, 2, 4), xaxt = "n", ...)
    dotchart(y, groups = groups, color = c(4, 2)[groups], 
             font = 2, pch = 19, gcolor = c(4, 2), xlab = "Participants' Gain Scores",
             pt.cex = ifelse(n1 <= 20 || n2 <= 20, 1.5, .8), labels = NA, main = NA)
  }
  par(xaxt = "s") ; axis(1, font = 2)
  
  gpos = rev(cumsum(rev(tapply(groups, groups, length)) + 2) - 1)
  
  u = par("usr")  
  
  segments(c(mean.g2, mean.g1), c(u[3], u[4]), c(mean.g2, mean.g1), rep(gpos[[2]], 2), lty = 2,
           col = c(2, 4))
  
  arrows(mean.g2, gpos[[2]], mean.g1, gpos[[2]], code = 3, length = .08, col = "darkgreen")
  
  mean.diff = mean.g1 - mean.g2
  
  text((mean.g1+mean.g2)/2, gpos[[2]], bquote(bold("Mean diff." == .(decimal((mean.diff), 2)))), font = 2, pos = 3, col = "green4", cex = 1.15 )
  
  legend("topright", legend = bquote(bold("Cohen's"~ bolditalic(d) == .(decimal(Cohend, 2)) )), bty = "n", text.col = "red4", cex = 1.15, bg = NA)
  
  if(descriptives) {
    
    legend("topleft", legend = bquote(bold(Mean == .(decimal(mean.g1, 2)))), text.col = 4, bty = "n", bg = NA)
    
    legend("topleft", legend = bquote(bold(sd == .(decimal(sd.g1, 2)))), text.col = 4, bty = "n", bg = NA,
           inset = .03, adj =  c(.2, 0.5) )
    
    legend("bottomleft", legend = bquote(bold(Mean == .(decimal(mean.g2, 2)))), text.col = 2, bty = "n", bg = NA, 
           inset = .03, adj = .1)
    legend("bottomleft", legend = bquote(bold(sd == .(decimal(sd.g2, 2)))), text.col = 2, bty = "n", bg = NA,
           adj =  c(-.1, 0.5))
  }
  m = matrix(c("R", "R", "O1", "O3", "T", "", "O2", "O4", "->", "->", "O2 - O1", "O4 - O3", "->", "->", "GainT", "GainC"), nrow = 2)
  dimnames(m) = list("PRE-POST-CONTROL DESIGN:" = c("", ""), c(rep("", 8)))
  m <- noquote(m)
  
  test <- data.frame(test[1:3], row.names = "result:")
  colnames(test) <- c("t.value", "df", "p.value")
  print(list(m, test), digits = digits)
}
# EXAMPLE OF USE:
pre.post.cont(n1 = 30, n2 = 30, effect.size = .5)
