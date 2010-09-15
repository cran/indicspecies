plotcoverage<-function(speciescomb, by=0.05,...) {
   Atseq <- seq(0,1,by=by)
   A <- speciescomb$A
   ci<-FALSE
   if(is.data.frame(A)) {
   	  LA<-A$lowerCI
   	  UA<-A$upperCI
   	  A<-A$stat
      covLA<-covUA<-numeric(length(Atseq))
      ci<-TRUE
   } 
   covA<-numeric(length(Atseq))
   for(i in 1:length(Atseq)) {
   	  sel = A>Atseq[i]
   	  sel[is.na(sel)]<-TRUE
   	  if(sum(sel)>0) covA[i] = coverage(speciescomb, selection = sel)
   	  else covA[i] = 0
   }
   plot(Atseq,covA*100, type="l", axes=FALSE, xlab=expression(A[t]), ylab="Coverage (%)",...)	
   if(ci) {
	   for(i in 1:length(Atseq)) {
    	  sel = LA>Atseq[i]
   	      sel[is.na(sel)]<-TRUE
   		  if(sum(sel)>0) covLA[i] = coverage(speciescomb, selection = sel)
   	  	  else covLA[i] = 0
    	  sel = UA>Atseq[i]
   	      sel[is.na(sel)]<-TRUE
   		  if(sum(sel)>0) covUA[i] = coverage(speciescomb, selection = sel)
   	  	  else covUA[i] = 0
   	    }
   	    lines(Atseq, covLA*100, lty=2, col="gray")
   	    lines(Atseq, covUA*100, lty=2, col="gray")
   }
   axis(1, pos=0)
   axis(2, pos=0)
}
