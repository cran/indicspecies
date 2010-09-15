#Finds the combination of clusters which is most significantly associated to each of the species patterns
multipatt <- function (x, cluster, func = "IndVal.g", duleg = FALSE, restcomb=NULL, nperm = 999,                                                                  
  torus = FALSE, grid.size, print.perm = FALSE)                                                                                              
{
	                                                                                                                              

# Matrix of possible cluster combinations
cl.comb <- function(clnames) {
	  k <- length(clnames)
    ep <- diag(1,k,k)
    names.ep <- clnames
    for(j in 2:k) {
      nco <- choose(k,j)
      co <- combn(k,j)
      epn <- matrix(0,ncol=nco,nrow=k)
      for(i in 1:ncol(co)) {
	  epn[co[,i],i] <- 1
	  names.ep <- c(names.ep, paste(clnames[co[,i]], collapse = "+"))
      }
      ep <- cbind(ep,epn)
    }
    colnames(ep) <- names.ep
    return(ep)
}

# Correlation measures for combinations
rcomb <- function(x, cluster, comb, k, mode="group", duleg=FALSE, restcomb=NULL) {
      nsps = ncol(x)
      N = length(cluster)	
      ni = diag(t(comb) %*% comb)[1:k]
      tx <- t(x)
      aisp = (tx %*% comb)[,1:k]
      lisp = (tx^2 %*% comb)[,1:k]

  if(mode=="site") {
      lspK = rowSums(lisp)
      aspK = rowSums(aisp)		
      aspC = (tx %*% comb)
      nC = diag(t(comb) %*% comb)
  } else if(mode=="group") {
      aispni=sweep(aisp,2,ni,"/")
      lispni=sweep(lisp,2,ni,"/")
      #Corrected sum and length of species vectors
      lspK = (N/k) * rowSums(lispni)
      aspK = (N/k) * rowSums(aispni)

      if(duleg) {
	  aspC = (N/k)*aispni
	  nC = rep(N/k,k)
      } else {
	  #Corrected sum of species values in combinations
	  aspC = matrix(0,nrow=nsps,ncol=ncol(comb))
	  #Corrected size of cluster combinations
	  nC = vector(mode="numeric",length=ncol(comb))
	  #Level 1
	  aspC[,1:k] = aispni[,1:k]
	  nC[1:k] = 1
	  #Remaining levels
	  cnt = k+1	   
	  for(level in 2:k) {
	    co = combn(1:k,level)
	    for(j in 1:ncol(co)) {
		aspC[,cnt] = rowSums(aispni[,co[,j]])
		nC[cnt] = length(co[,j])
		cnt= cnt+1
	    }
	  }
	  aspC = (N/k)*aspC
	  nC = (N/k)*nC
      }
  }

  #Compute index
  num = N*aspC - aspK%o%nC
  den = sqrt(((N*lspK)-aspK^2)%o%(N*nC-(nC^2)))
  str=num/den
  if(!duleg) str <- str[,-ncol(str)] # remove all sites as combination for correlation indices
  if(!duleg && !is.null(restcomb)) str <- str[,restcomb]
#  colnames(str) <- colnames(comb)[1:ncol(str)]
  return(str)
}

# IndVal for combinations
indvalcomb <- function(x, cluster, comb, k, mode = "group", duleg = FALSE, restcomb=NULL) {
  tx <- t(x)
  aisp = tx %*% comb
  dx <- dim(tx)
  nisp <- matrix(as.logical(tx),nrow=dx[1],ncol=dx[2]) %*% comb
#    nisp = t(ifelse(x > 0, 1, 0)) %*% comb
  ni = diag(t(comb) %*% comb)
  nispni = sweep(nisp, 2, ni, "/")   
  if (mode == "site") A = sweep(aisp, 1, colSums(x), "/")  
  else {
    aispni = sweep(aisp[, 1:k], 2, ni[1:k], "/")
    asp = rowSums(aispni[, 1:k])
    if(duleg) A = sweep(aispni, 1, asp, "/") 
    else {
    s = aispni #matrix(0, nrow = nrow(aispni), ncol = ncol(comb))
    for(j in 2:(k)) {
      co <- combn(k,j)
      sn <- apply(co, 2, function(x) rowSums(aispni[,x]))
      s <- cbind(s, sn)
    }
	  A = sweep(s, 1, asp, "/")
  	} 
  }
  iv = sqrt(A * nispni)
  colnames(iv) <- colnames(comb)
  if(!duleg && !is.null(restcomb)) iv = iv[,restcomb]
  return(iv)
}


sampletorus<-function(grid.size) {
  d1= grid.size[1]
  d2= grid.size[2]
  m1 = as.integer(runif(1,0,d1))
  m2 = as.integer(runif(1,0,d2))
  nsites= d1*d2
  pInd = 1:nsites
  for(j in 1:nsites) {
	  yc = (j-1)%%d2+1
	  xc = (j-yc)/d2
	  r = ((xc+m1)%%d1)*d2 + (yc-1+m2)%%d2 + 1
	  pInd[j]=r
  }
  return(pInd)
}
	
  vegnames <- names(x)
  x <- as.matrix(x)                                                                                                                                
  nsps = ncol(x)
  clnames = levels(as.factor(cluster))
  k = length(clnames)

  # Check parameters
   func= match.arg(func, c("r","r.g","IndVal.g","IndVal"))
	if(torus && is.null(grid.size)) stop("Please, supply 'grid.size' dimensions if you want to use torus translation.")
	if(k<2) stop("At least two clusters are required.")
  if(sum(is.na(cluster))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(sum(is.na(x))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(!is.null(restcomb)) restcomb = as.integer(restcomb)

  # creates combinations from clusters
  combin <- cl.comb(clnames)	# possible combinations (can also be used for permutations)
  if (duleg) combin <- combin[,1:k]
  
  #Builds the plot membership matrix corresponding to combinations
  comb <- combin[cluster,]

  # Computes association strength for each group
  str <- switch(func,
    r =		rcomb(x, cluster, comb, k, mode = "site", duleg = duleg, restcomb = restcomb),
    r.g = 	rcomb(x, cluster, comb, k, mode = "group", duleg = duleg, restcomb = restcomb),
    IndVal = 	indvalcomb(x, cluster, comb, k, mode = "site", duleg, restcomb = restcomb),
    IndVal.g = 	indvalcomb(x, cluster, comb, k, mode = "group", duleg, restcomb = restcomb),
    )

  # Maximum association strength
  maxstr = apply(str,1,max) 
  wmax <- max.col(str)
  #prepares matrix of results
  if(!duleg && !is.null(restcomb))  m <- as.data.frame(t(combin[,restcomb][,wmax]))
  else  m <- as.data.frame(t(combin[,wmax]))
  dimnames(m) <- list(vegnames, sapply(clnames, function(x) paste("s", x, sep='.')))
  m$index <- wmax
  m$stat <- apply(str,1,max)

  #Perform permutations and compute p-values
  pv <- 1
  for (p in 1:nperm) {
      if (!torus) pInd <- sample(1:length(cluster)) else pInd <- sampletorus(grid.size)
      tmpcls = cluster[pInd]
		 combp = combin[tmpcls,]
      tmpstr <- switch(func,
	r   = rcomb(x, cluster = tmpcls, combp, k, mode = "site", duleg = duleg, restcomb = restcomb),
	r.g = rcomb(x, cluster = tmpcls, combp, k, mode = "group", duleg = duleg, restcomb = restcomb),
	IndVal = indvalcomb(x, cluster = tmpcls, combp, k, mode = "site", duleg, restcomb = restcomb),
	IndVal.g= indvalcomb(x, cluster = tmpcls, combp, k, mode = "group", duleg, restcomb = restcomb)
      )
      tmpmaxstr <- vector(length=nrow(tmpstr))
      for(i in 1:nrow(tmpstr)) tmpmaxstr[i] <- max(tmpstr[i,])	# apply is more slowly in this case
      pv = pv + (tmpmaxstr >= m$stat)
  }

  m$p.value <- pv/(1 + nperm)
  #Put NA for the p-value of species whose maximum associated combination is the set of all combinations
  m$p.value[m$index == (2^k-1)] <- NA

  if(!duleg && !is.null(restcomb))  comb<-comb[,restcomb]
  
  a = list(call=match.call(), func = func, cluster = cluster, comb = comb, str = str, sign = m)
  class(a) = "multipatt"
  return(a)
}


