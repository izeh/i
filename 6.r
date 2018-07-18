source("https://raw.githubusercontent.com/rnorouzian/i/master/i.r")
C = beta.id(.6, .8, .6)
prior = function(x) dbeta(x, C$a, C$b) 
likelihood = function(x) dbinom(55, 100, x)
posterior = function(x) prior(x)*likelihood(x)
mode = optimize(posterior, 0:1, maximum = TRUE)[[1]]

par(mfcol = c(1, 3), mar = c(5, .4, 15, .1), xpd = NA, yaxt = "n", bty = "n", yaxs = "i", lwd = 2, las = 1, font.lab = 2, font = 2)

Prior = curve(prior, n = 1e4, xlab = "Proportion of (B)", ylab = NA, xaxt = "n")$y
mtext("\u00D7", 4, line = 1, cex = 3, font = 1)
text(.7, max(Prior), "Prior", pos = 3)
axis(1, at = x <- axTicks(1), lab = paste0(x*1e2, "%"))

Likelihood = curve(likelihood, n = 1e4, xlab = "Proportion of (B)", ylab = NA, xaxt = "n")$y
mtext("=", 4, line = 1, cex = 3, font = 1)
text(.55, max(Likelihood), "Likelihood", pos = 3)
axis(1, at = x <- c(0, .2, .4, .7, .85, 1), labels = paste0(x*1e2, "%"))
axis(1, at = .55, labels = paste0(55, "%"), col = 2, col.axis = 2)

Posterior = curve(posterior, n = 1e4, xlab = "Proportion of (B)", ylab = NA, xaxt = "n")$y
text(.58, max(Posterior), "Posterior", pos = 3)
axis(1, at = x <- c(0, .2, .4, .7, .85, 1), labels = paste0(x*1e2, "%"))
axis(1, at = mode, labels = paste0(round(mode, 4)*1e2, "%"), col = 2, col.axis = 2)
