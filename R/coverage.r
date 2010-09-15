	#Determines the coverage of a set of indicators
	coverage <- function (speciescomb, selection=NULL) {
		if(is.null(selection)) selection = rep(TRUE, nrow(speciescomb$C))
		if(length(dim(speciescomb$C))==2) c = speciescomb$C[selection,]
		else c = speciescomb$c[selection]
		group.vec = speciescomb$group.vec
		if(length(dim(speciescomb$XC))==2) xc = speciescomb$XC[, selection]
		else xc = speciescomb$XC[selection]
		if(sum(selection)>1) {
			ccx = rep(FALSE, nrow(xc))    
			for(rc in 1:nrow(c)) {
  	  			ccx = ccx | xc[,rc]>0  		
  			}
  			return(sum(ccx & group.vec) / sum(group.vec))
  		} else {
  			return(sum(xc>0 & group.vec) / sum(group.vec))
  		}
  	}