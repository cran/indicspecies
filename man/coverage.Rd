\encoding{UTF-8}
\name{coverage}
\alias{coverage}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determines the coverage of a set of indicators}
\description{
 Calculates the proportion of sites of the target site group where one or another indicator (species combination) is found.
}
\usage{
coverage(speciescomb, selection=NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{speciescomb}{An object of class '\code{\link{speciescomb}}'.}
  \item{selection}{A logical vector restricting the set of indicators used to calculate the coverage.}
}
\value{
  The proportion of sites of the target site group where one or another indicator (species combination) is found.
}
\references{
De \enc{Cáceres}{Caceres}, M., Legendre, P., Wiser, S.K. & Brotons, L (in preparation). Using species combinations in indicator analyses.

De \enc{Cáceres}{Caceres}, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.

}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Biodiversity and Landscape Ecology Laboratory, Forest Science Center of Catalonia
}
\seealso{\code{\link{speciescomb}}, \code{\link{pruneindicators}}}
\examples{
library(stats)

data(wetland) ## Loads species data

## Creates three clusters using kmeans
wetkm = kmeans(wetland, centers=3) 

## Determine sensitivity of individual species
B=strassoc(wetland, cluster=wetkm$cluster,func="B") 

## Select species with more than 20% of sensitivity for the first group
sel=which(B[,1]>0.2) 

## Run indicator analysis with species combinations for the first group
sc=speciescomb(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)

## Determine the coverage of the selected set of indicators
coverage(sc)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.