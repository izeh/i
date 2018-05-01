source("https://raw.githubusercontent.com/rnorouzian/i/master/i.r")
C = beta.id(.6, .8, .6) ; options(warn = -1)
curve(dbeta(x, C[1,], C[2,]), n = 1e4, axes = FALSE, lwd = 2, yaxs = "i", xpd = TRUE, 
      xlab = "Proportion of preference for (B)", ylab = NA, font.lab = 2)

axis(1, at = axTicks(1), labels = paste0(axTicks(1)*1e2, "%") )
axis(1, at = axTicks(1)[4:5], labels = paste0(axTicks(1)[4:5]*1e2, "%"), font = 2)
x = seq(.6, .8, l = 10)
y = dbeta(x, C[1,], C[2,])
arrows(x, 0, x, y, code = 2, length = .13, angle = 10, col = 4)
x = c(seq(0, .577, l = 28), seq(.822, 1, l = 10))
y = dbeta(x, C[1,], C[2,])
arrows(x, 0, x, y, code = 1, length = .12, angle = 10, col = 2)
