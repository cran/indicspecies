summary.multipatt <- function (object, alpha = 0.05, minstat = 0, indvalcomp=FALSE,...) {
    x <- object
    ncomb = ncol(x$str)
    ncolsign = ncol(x$sign)
    nsps = nrow(x$sign)
    cat("\n Multilevel pattern analysis")
    cat("\n ---------------------------\n")
    cat("\n Association function:", x$func)
    cat("\n Significance level (alpha):", alpha)
    cat("\n Minimum statistic value (minstat):", minstat, "\n")
    cat("\n Total number of species:", nsps)
    a = x$sign[!is.na(x$sign$p.value) & x$sign$p.value <= alpha & 
	x$sign$stat > minstat, ]
    cat("\n Selected number of species:", nrow(a), "\n")
    cols = (ncolsign - 1):ncolsign
    
	#Collate indval components if required and possible
	if(indvalcomp && !is.null(x$B) && !is.null(x$A)) {
		As = numeric(nrow(x$sign))
		Bs = numeric(nrow(x$sign))
		for(i in 1:nrow(x$sign)) {
			As[i] = x$A[i,x$sign$index[i]]
			Bs[i] = x$B[i,x$sign$index[i]]
		}
		y = cbind(x$sign, As,Bs)
		cols = c(ncol(y)-1, ncol(y),cols)
		names(y) = c(names(x$sign),"A","B")
	}
	else y = x$sign
  
  #Show indicators for all combinations  
  for (i in 1:ncomb) {
	sum = 1
	for (k in 1:(ncolsign - 4)) {
	    if (i == sum) 
		cat("\n Number of species associated to", k, if(k==1) "group:" else "groups:",
		 sum(rowSums(a[, 1:(ncolsign-3)]) == k), "\n")
	    sum = sum + choose(ncolsign-3, k)
	}
	m = y[x$sign$index == i & !is.na(x$sign$p.value) & 
	    x$sign$p.value <= alpha & x$sign$stat > minstat, ]
	if (nrow(m) > 0) {
	    cat("\n Group", dimnames(x$comb)[[2]][i], " #sps. ", nrow(m), "\n")
	    m = m[order(m$stat, decreasing = TRUE), cols]
	    printCoefmat(m, signif.stars = TRUE, signif.legend = FALSE, 
		digits = 4, P.values = TRUE, has.Pvalue = TRUE)
	}
    }
    Signif <- symnum(x$sign$p.value, corr = FALSE, na = FALSE, 
	cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
    cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
}
