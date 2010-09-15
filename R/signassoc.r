`signassoc` <-
function(X, U=NULL, cluster=NULL, mode = 1, nperm=999, alternative="greater", torus=FALSE, grid.size, print.perm=FALSE) {
	
vector.to.partition <- function(v, clnames) {
    m <- t(sapply(v,function(x) as.numeric(x==clnames)))
    dimnames(m) = list(1:length(v),clnames)
    return(m)                                                                                                                              
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

  nsps = ncol(X)
  nsites = nrow(X)
  
   mode= match.arg(as.character(mode), c("0","1"))
   alternative= match.arg(as.character(alternative), c("greater","less","two.sided"))
	if(torus && is.null(grid.size)) stop("Please, suply 'grid.size' dimensions if you want to use torus translation.")
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  
  if(is.null(U)) U = vector.to.partition(cluster, levels(factor(cluster))) 
  if(sum(is.na(U))>0) stop("Cannot deal with NA values. Remove and run again.")
  
  ngroups = ncol(U)	
  cdm = matrix(1,nrow=nsps,ncol=ngroups)
  ddm = matrix(1,nrow=nsps,ncol=ngroups)
  
  X = as.matrix(X)
  U = as.matrix(U)
  if(mode==0) {
	  dm = t(X)%*%U
  } else {
	  aisp = t(X)%*%U
  	  ni = diag(t(U)%*%U)
  	  aispni=sweep(aisp,2,ni,"/")
  	  aispni[is.na(aispni)]=0 # check for division by zero
  	  s = apply(aispni,1,"sum")
  	  dm = sweep(aispni,1,s,"/")
  	  dm[is.na(dm)]=0 # check for division by zero
  }
  	
  a <- system.time({
  for(p in 1:nperm) {
      if(p%%100==0 & print.perm) cat("perm", p,"\n")
		if(!torus) {
		   pInd=sample(1:nsites)
		} else {
		   pInd=sampletorus(grid.size)
		}
		pX = as.matrix(X[pInd,])
  		if(mode==0) {
	  		dmp = t(pX)%*%U
  		} else {
	 		aisp = t(pX)%*%U
		  	ni = diag(t(U)%*%U)
  	  		aispni=sweep(aisp,2,ni,"/")
	    	aispni[is.na(aispni)]=0 # check for division by zero
  	  		s = apply(aispni,1,"sum")
  	  		dmp = sweep(aispni,1,s,"/")
	   	   dmp[is.na(dmp)]=0 # check for division by zero
  		}			
		if(alternative=="less") {
		   cdm = cdm + as.numeric(dmp<=dm)
		} else if(alternative=="greater") {
		   cdm = cdm + as.numeric(dmp>=dm)
		} else if(alternative=="two.sided") {
		   cdm = cdm + as.numeric(dmp>=dm)
		   ddm = ddm + as.numeric(dmp<=dm)
		}
	}
	if(alternative!="two.sided") {
		cdm=cdm/(nperm+1)
	} else {
		cdm=pmin(matrix(1,nrow=nsps,ncol=ngroups),(2*pmin(cdm,ddm))/(nperm+1))
	}
	cdm = as.data.frame(cdm)
	row.names(cdm)=names(X)
   names(cdm)=names(as.data.frame(U))
	})
  a[3] <- sprintf("%2f",a[3])
  #cat("Time to compute p-values =",a[3]," sec",'\n')


   psidak = vector(mode="numeric", length=nsps)
   best = vector(mode="numeric", length=nsps)
   for(i in 1:nsps) {
   	  best[i] = which(cdm[i,]==min(cdm[i,]))[1]
	  psidak[i] = (1-(1-min(cdm[i,]))^ngroups)
   	}
	cdm = data.frame(cbind(cdm,best,psidak))
	names(cdm) = c(levels(as.factor(cluster)),"best","psidak")
	return(cdm)      
 }

