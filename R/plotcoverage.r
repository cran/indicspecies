plotcoverage<-function(x, by=0.05, type="stat", max.order=NULL, add=FALSE, ...) {
   match.arg(type,c("lowerCI","upperCI","stat"))
   Atseq <- seq(0,1,by=by)
   A <- x$A
   if(is.data.frame(A)) {
      if(type=="lowerCI") A<-A$lowerCI
   	  else if(type=="upperCI") A<-A$upperCI
   	  else A<-A$stat
      covLA<-covUA<-numeric(length(Atseq))
   } 
   num.order <- rowSums(x$C)
   sel2 <- rep(TRUE,length(num.order))
   if(!is.null(max.order)) sel2[num.order>max.order] <- FALSE
   covA<-numeric(length(Atseq))
   for(i in 1:length(Atseq)) {
   	  sel = A>Atseq[i]
   	  sel[is.na(sel)]<-TRUE
   	  if(sum(sel & sel2)>0) covA[i] = coverage(x, selection = sel & sel2)
   	  else covA[i] = 0
   }
   if(!add) {
     plot(Atseq,covA*100, type="s", axes=FALSE, xlab=expression(A[t]), ylab="Coverage (%)", ...)	
     axis(1, pos=0)
     axis(2, pos=0)     
   } else lines(Atseq,covA*100, type="s", ...)
}
