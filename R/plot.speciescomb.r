plot.speciescomb<-function(x, type="sqrtIV", maxline=TRUE, ...) {
   A = x$A
   B = x$B
   order = rowSums(x$C)
   if(is.data.frame(A)) {
	   if(type=="IV") val = A[,1] * B[,1]
   	   else if(type=="sqrtIV") val = sqrt(A[,1] * B[,1])
   	   else if(type=="A") val = A[,1]	
       else if(type=="B") val = B[,1]	
       else if(type=="LA") val = A[,2]	
       else if(type=="UA") val = A[,3]	
       else if(type=="LB") val = B[,2]	
       else if(type=="UB") val = B[,3]	
   } else {
	   if(type=="IV") val = A * B
   	   else if(type=="sqrtIV") val = sqrt(A * B)
   	   else if(type=="A") val = A	
       else if(type=="B") val = B	   	
   }
   plot(order,val, type="n", axes=FALSE, xlab="Order", ylab=type,...)	
   points(order,val, pch=1, cex=0.5)	
   axis(1)
   axis(2)
   if(maxline) lines(1:max(order),tapply(val,order,max), col="gray")
}
