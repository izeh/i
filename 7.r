source("https://raw.githubusercontent.com/rnorouzian/i/master/i.r")
C = beta.id(.6, .8, .6)
prior = function(x) dbeta(x, C$a, C$b) 
likelihood = function(x) dbinom(55, 100, x)
posterior = function(x) prior(x)*likelihood(x)

CI = HDI(posterior)

Posterior = curve(posterior, n = 1e4, axes = FALSE, lwd = 2, yaxs = "i",
                  xlab = "Proportion of preference for (B)", ylab = NA, font.lab = 2)$y

axis(1, at = x <- c(0, .2, .4, .7, .85, 1), labels = paste0(x*1e2, "%"), font = 2)
mode = optimize(posterior, 0:1, maximum = TRUE)[[1]]
axis(1, at = mode, labels = paste0(round(mode, 4)*1e2, "%"), font = 2, col = 2, col.axis = 2)
segments(mode, 0, mode, max(Posterior), lty = 3, col = "darkgreen")

segments(CI[1], 0, CI[2], 0, lend = 1, lwd = 4, col = "magenta", xpd = NA)
points(mode, 0, pch = 21, bg = "cyan", cex = 1.5, col = "magenta", xpd = NA)
text(CI, 0, paste0(round(CI, 4)*1e2, "%"), pos = 3, font = 2, col = "magenta")
