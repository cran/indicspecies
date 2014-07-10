nichearea <-
function (P, D = NULL, axes=c(1,2)) {
    if (!inherits(P, "data.frame")) stop("Non convenient dataframe for species resource use")
    if (!is.null(D)) {
        if (!inherits(D, "dist")) stop("Object of class 'dist' expected for distance")
        D <- as.matrix(D)
        if (ncol(P) != nrow(D)) stop("The number of columns in P must be equal to the number of items in D")
        D <- as.dist(D)
    }
    #If no distance matrix is provided, the distance between resources is assumed to be maximum
    if (is.null(D)) D <- as.dist((matrix(1, ncol(P), ncol(P)) - diag(rep(1, ncol(P)))))
	 cmd = cmdscale(D,eig=TRUE,k= ncol(P)-1)
    X = cmd$points
    
    V <- as.data.frame(rep(0, nrow(P)))
    names(V) <- "Area"
    rownames(V) <- row.names(P)
    for (i in 1:nrow(P)) {
    	  pi = as.numeric(P[i,])
        if (is.na(sum(pi))) V[i, ] <- NA
        else if (sum(pi) < 1e-16) V[i, ] <- 0
        else {
   	   		if(sum(pi>0)==1) {
   	   			V[i,]=0
   	   		}else {
   	   			a =X[pi>0,axes]
	   	   		V[i,]=area.poly(as(a[chull(a),],"gpc.poly"))
   	   		}        	
        	}
    }
    return(V)
}

