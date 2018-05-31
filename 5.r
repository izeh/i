Likelihood = curve(dbinom(55, 100, x), n = 1e4, lwd = 2, axes = FALSE, yaxs = "i",
      xlab = "Proportion of preference for (B)", ylab = NA, font.lab = 2)$y

axis(1, at = x <- c(0, .2, .4, .7, .85, 1), labels = paste0(x*1e2, "%"), font = 2)
axis(1, at = .55, labels = paste0(55, "%"), font = 2, col = 2, col.axis = 2)
segments(.55, 0, .55, dbinom(55, 100, .55), lty = 3, col = 2)
text(.55, max(Likelihood)/2, "Obtained Proportion of (B)", pos = 3, srt = 90, font = 2, col = "red4", cex = .9)
