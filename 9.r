source("https://raw.githubusercontent.com/rnorouzian/i/master/i.r")

d.ci.sim <- function(d = .5, n1 = 30, n2 = NA, conf.level = .95, n.sim = 5, ylabel = TRUE){

fun <- function(){  
  ds <- rcohen(1, dbase = d, n1 = n1, n2 = n2)
  c(CI = as.numeric(d.ci(d = ds, n1 = n1, n2 = n2, conf.level = conf.level)[,2:3]), ds = ds)
}

sim <- t(replicate(n.sim, fun()))
capture <- sim[ ,1] <= d & d <= sim[ ,2]

par(mgp = c(2, .2, 0), tck = -.015)    
plot(sim[, 1:2], rep(1:n.sim, 2), ty = "n", ylab = NA, yaxt = "n", xlab = "Effect Size (Cohen's d)", font.lab = 2)
axis(1, at = d, col.axis = 2, col = 2, font = 2)
abline(h = 1:n.sim, col = 8, lty = 3)
if(ylabel) axis(2, at = 1:n.sim, labels = paste0("Repeat ", rev(1:n.sim)), font = 2, las = 1, cex.axis = .6, tck = -.006)
abline(v = d, lty = 2, col = 2) 
segments(sim[ ,1], 1:n.sim, sim[ ,2], 1:n.sim, lend = 1, col = ifelse(capture, 1, 2))
points(sim[, 3], 1:n.sim, pch = 19, col = ifelse(capture, 1, 2), cex = ifelse(n.sim > 50, .6, .65))

cat(paste0("\t", "Coverage = ", mean(capture)*1e2, "%")) 
}
# Example of use:
d.ci.sim(d = .5, n1 = 30, n2 = NA, n.sim = 20)
