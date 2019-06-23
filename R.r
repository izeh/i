## helper function:
sdif <- function(n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, t = NA, F1 = NA, sdp = NA){
  
  ifelse(!is.na(r) & !is.na(sdpre) & !is.na(sdpos), sqrt(sdpre^2+sdpos^2-2*r*sdpre*sdpos),
         ifelse(!is.na(n) & is.na(r) & !is.na(t) & !is.na(mpre) & !is.na(mpos), sqrt((n*(mpos - mpre)^2)/(ifelse(is.na(F1) & !is.na(t), t^2, ifelse(!is.na(F1) & is.na(t), F1, NA)))), 
                ifelse(!is.na(r) & !is.na(sdp), sqrt(2*sdp^2*(1-r)), NA)))
  }

## Actual function:
d.prepos <- function(study.name = NA, group.name = NA, n = NA, mpre = NA, mpos = NA, sdpre = NA, sdpos = NA, r = NA, autoreg = FALSE, t = NA, sdif = NA, sdp = NA, F1 = NA, df2 = NA, post, control, ...) 
{
  
  if(missing(control) || missing(post)) stop("'post' or/and 'control' missing.", call. = FALSE)  
  
  d <- NA
  mdif <- ifelse(!is.na(mpre) & !is.na(mpre), mpos - mpre, NA)
  sdif <- ifelse(is.na(sdif), sdif(sdpre = sdpre, sdpos = sdpos, t = t, r = r, n = n, mpos = mpos, mpre = mpre, F1 = F1, sdp = sdp), sdif)
  ifelse(!is.na(mdif) & is.na(d) & !is.na(sdif), d <- mdif/sdif, d <- NA) 
  
  out <- data.frame(d = d, n = n, sdif = sdif, post, control, ...)
  
  if(all(is.na(out$d))) stop("\ninsufficient info. to calculate effect size(s).", call. = FALSE)
  
  return(out)
}
########################################## YOUR SUGGESTED SOLUTIONS: ###########################################################
library(purrr)
D <- read.csv("https://raw.githubusercontent.com/izeh/i/master/k.csv", h = T)
L <- split(D, D$study.name) ; L[[1]] <- NULL


dot.names <- names(L[[1]])[!names(L[[1]]) %in% formalArgs(d.prepos)]

args <- lapply(L, function(x) unclass(x[c(head(formalArgs(d.prepos), -1), dot.names)]))

g <- simplify2array(lapply(names(args[[1]]), function(i) 
  lapply(args, function(j) j[i])))

do.call(Map, c(f = d.prepos, unname(split(g, col(g)))))        ## Fails <><><><><><><><><>

unname(do.call(Map, c(f = d.prepos, purrr::transpose(args))))  ## Works <><><><><><><><><>
