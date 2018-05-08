biling.survey <- function(N = 40, pBi = .5, ...){
  
  Responses = sample(0:1, size = N, prob = c(1 - pBi, pBi), replace = TRUE)  # Generate N responses from parents
  # (B = 1, M = 0)
  Prop = cumsum(Responses) / 1:N   # Compute proportion of B as each parent responds
  
  par(las = 1, tck = -.02, font.lab = 2, ...)  
  plot.ts(Prop, type = "o", ylim = 0:1, yaxt = "n", pch = 21, bg = 3, xlab = "Parent ID number", 
          ylab = "Proportion of (B)", xaxt = "n")
  
  CI = binom.test(cumsum(Responses)[N], N)[[4]]
  arrows(N, CI[1], N, CI[2], code = 3, angle = 90, length = .08, lend = 1)
  points(N, Prop[N], pch = 21, bg = "cyan")
  
  axis(1, at = axTicks(1), lab = c(1, axTicks(1)[-1]))  
  axis(2, at = seq(0, 1L, len = 6), lab = paste0(seq(0, 1e2, len = 6), "%"))
  
  ResponseSeq = paste(c("M", "B")[Responses[1L:1e1] + 1L], collapse = "")
  
  Display = paste0("Response Sequence = ", ResponseSeq, ". . .")
  CI = paste0("95% CI: [", round(CI[1], 4)*1e2, "%", ", ", round(CI[2], 4)*1e2, "%", "]")
  
  text(N, c(1, .95, .9), c(Display, paste0("Proportion of B in ", N, " responses = ", Prop[N]*1e2, "%"), CI), adj = c(1, .5), col = "red4", font = 2, cex = .8)  
}
#Example of use:
biling.survey(N = 100, pBi = .75)
