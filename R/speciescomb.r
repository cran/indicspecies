#For any combination of species within the given data calculate the probability that the plot observation
#belongs to the given type once the species in the combination have been found altogether.
speciescomb <- function (X, cluster, group, func="IndVal", max.order = 5, At = 0, Bt=0, sqrtIVt =0, nboot=0, alpha=0.05, XC = TRUE, verbose=FALSE) {
	                 
  func <- match.arg(func, c("IndVal", "IndVal.g"))                                                                                                             
  if(sum(is.na(cluster))>0) stop("Cannot deal with NA values. Remove and run again.")
  if(sum(is.na(X))>0) stop("Cannot deal with NA values. Remove and run again.")
  nsites = nrow(X)

  cluster= as.factor(cluster)
  group <- as.character(group)
  group <- match.arg(group, levels(cluster))                                                                                                             
  if(verbose) cat(paste("Target site group: ",group,"\n", sep=""))
  group.vec <- cluster==group
  group.vec[is.na(group.vec)] = FALSE

  #Get species names
  spplist = names(X)
  if(verbose) cat(paste("Number of candidate species: ",length(spplist),"\n", sep=""))
  if(length(spplist)==1) stop("At least two species are necessary.")
  
  #Select rows that contain the species or the group
  ng = sum(group.vec)
  if(verbose) cat(paste("Number of sites:",nsites,"\n"))
  if(verbose) cat(paste("Size of the site group:",ng,"\n"))
  
  #Create structures to store data
  Cvalid<-matrix(0,nrow=0,ncol=length(spplist))
  Astat = numeric(0)
  Bstat = numeric(0)
  sqrtIVstat = numeric(0)
  k = length(spplist)
  totco = 0 
  for(j in 1:min(max.order, k)) {
      co <- combn(k,j)
      totco = totco + ncol(co)
	  if(verbose) cat(paste("Evaluating ", ncol(co) ," combinations of ",j," species",sep=""))
	  A=rep(NA,ncol(co))
  	  B=rep(NA,ncol(co))
	  for(r in 1:ncol(co)) {
  		  if(ncol(co)>100) if(r%%round(ncol(co)/10)==0 && verbose) cat(".")
  		  if(j==1) sc.ab<-X[,co[,r]]
  		  else sc.ab<-apply(X[,co[,r]],1,min)
  		  if(sum(sc.ab)>0) {
  			scg = sc.ab[group.vec]
  			if(func=="IndVal.g") {
  				mg = (sum(scg)/ng)
		  		A[r] = mg/sum(tapply(sc.ab,cluster, "mean"))
  			} else {
  				A[r] = sum(scg)/sum(sc.ab)
  			}
  			B[r] = sum(scg>0)/ng
    	  } 
  	  }
  	  sqrtIV = sqrt(A*B)
	  #Remove non-valid combinations
 	  sel = A>=At & B>=Bt & sqrtIV>=sqrtIVt
	  #Remove combinations that do not have any co-occurrence (NAs)
 	  sel[is.na(sel)] = FALSE
	  if(verbose) cat(paste(" - ", sum(sel), " valid combinations.\n", sep=""))
	  if(sum(sel)>0) {
   	   	epn <- matrix(0,nrow=sum(sel),ncol=k)
   	   	indices = which(sel)
      	for(i in 1:length(indices)) epn[i,co[,indices[i]]] <- 1
	  	Astat = c(Astat,A[sel])
	  	Bstat = c(Bstat,B[sel])
	  	sqrtIVstat= c(sqrtIVstat,sqrtIV[sel])
	  	Cvalid = rbind(Cvalid,epn)
	  	rm(epn)
	  }
	  rm(co)
	  rm(A)
	  rm(B)
	  rm(sqrtIV)
  	  gc()
  }
  if(verbose) cat(paste("Number of combinations explored: ",totco,"\n", sep=""))
  if(verbose) cat(paste("Number of valid combinations: ",nrow(Cvalid),"\n", sep=""))
  if(nrow(Cvalid)==0) return()
  
  Cvalid = as.data.frame(Cvalid)
  names(Cvalid)<-spplist
  
  nc <- nrow(Cvalid)
 
  #Remove species that do not appear in any valid combination
  selSpp = colSums(Cvalid)>0
  if(verbose) cat(paste("Number of remaining species:",sum(selSpp),"\n"))
  Cvalid = Cvalid[,selSpp]
  X = X[,selSpp]
  nspp <- sum(selSpp)
  
  
  #Build abundance matrix for valid species combinations
  if(XC || nboot>0) {
  	if(verbose) {
  	 cat(paste("Building abundance matrix"))
  	}
  
  	if(nc>1) {
  		XC = data.frame(matrix(0, nrow=nrow(X), ncol=nc))
  		names(XC) =  row.names(Cvalid)
  	}
  	else XC = data.frame(rep(0,length(nrow(X))))
  	for(r in 1:nc) {
  		if(nc>100) if(r%%round(nrow(Cvalid)/10)==0 && verbose) cat(".")
  		if(nc>1 && nspp>1) {
  			if(sum(Cvalid[r,])==1) XC[,r]<-X[,Cvalid[r,]==1]
	  		else XC[,r]<-apply(X[, Cvalid[r,]==1],1,min)
	  	} else if(nc==1 && nspp>1) {
	  		if(sum(Cvalid)==1) XC<-X[,Cvalid==1] 
	  		else XC<-apply(X[, Cvalid==1],1,min)
	  	} else if(nc==1 && nspp==1) {
	  		XC<-X 
	  	}
  	}
  	if(verbose) cat(paste("\n"))
  } else {
  	XC = NULL
  }
  
  #Calculate bootstrap confidence intervals for sensitivity and ppp of valid combinations
  if(nboot>0) {
  	  if(nboot<100) nboot=100 #Minimum of 100 bootstrap replicates
	  if(verbose) {
  			cat(paste("Calculating bootstrap confidence intervals"))
  	  }
	  dmbA = matrix(NA,nrow=nboot,ncol=nc)
	  dmbB = matrix(NA,nrow=nboot,ncol=nc)
	  dmbIV = matrix(NA,nrow=nboot,ncol=nc)
	  for(b in 1:nboot) {
	  	  if(b%%round(nboot/10)==0 && verbose) cat(".")
		  bi = sample(nsites,replace=TRUE)
		  ngb = sum(group.vec[bi])
		  XCB = as.matrix(XC[bi,])
		  XCBg = as.matrix(XCB[group.vec[bi],])
	  	  if(func=="IndVal.g") {
	  	  	kk <- colSums(apply(XCB,MARGIN=2,FUN=tapply,cluster[bi],"mean"))
		  	dmbA[b,] = colMeans(XCBg)/kk
  		  } else {
  			dmbA[b,] = colSums(XCBg)/colSums(XCB)
  		  }
		  dmbB[b,] = colSums(as.matrix(ifelse(XCBg>0,1,0)))/ngb
		  dmbIV[b,] = sqrt(dmbA[b,]*dmbB[b,])
	  }
	  if(verbose) cat(paste("\n"))
	  dmlowerA = rep(0,nc)
	  dmupperA = rep(0,nc)
	  dmlowerB = rep(0,nc)
	  dmupperB = rep(0,nc)
	  dmlowerIV = rep(0,nc)
	  dmupperIV = rep(0,nc)
	  for(i in 1:nc) {	
			sdmb = sort(dmbA[,i])			
			dmlowerA[i]=sdmb[(alpha/2.0)*nboot]
			dmupperA[i]=sdmb[(1-(alpha/2.0))*nboot]
			sdmb = sort(dmbB[,i])
			dmlowerB[i]=sdmb[(alpha/2.0)*nboot]
			dmupperB[i]=sdmb[(1-(alpha/2.0))*nboot]
			sdmb = sort(dmbIV[,i])
			dmlowerIV[i]= sdmb[(alpha/2.0)*nboot]
			dmupperIV[i]= sdmb[(1-(alpha/2.0))*nboot]
	  }
	  sA = as.data.frame(cbind(Astat,dmlowerA,dmupperA))
  	  names(sA) = c("stat", "lowerCI", "upperCI")
  	  row.names(sA) = row.names(Cvalid)
	  sB = as.data.frame(cbind(Bstat,dmlowerB,dmupperB))
  	  names(sB) = c("stat", "lowerCI", "upperCI")
  	  row.names(sB) = row.names(Cvalid)
	  sIV = as.data.frame(cbind(sqrtIVstat,dmlowerIV,dmupperIV))
  	  names(sIV) = c("stat", "lowerCI", "upperCI")
  	  row.names(sIV) = row.names(Cvalid)
  } else{
  	  sA = Astat
  	  sB = Bstat
  	  sIV = sqrtIVstat
  }
  
  result = list(C=Cvalid, XC=XC, A=sA, B=sB, sqrtIV=sIV, group.vec =group.vec)
  class(result) = "speciescomb"
  return(result)
}


