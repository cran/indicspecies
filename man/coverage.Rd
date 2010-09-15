\encoding{UTF-8}
\name{coverage}
\alias{coverage}
\alias{plotcoverage}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coverage of a set of indicators}
\description{
 Function \code{coverage} calculates the proportion of sites of the target site group where one or another indicator (species combination) is found. Function \code{plotcoverage} plots the coverage against the positive predictive value threshold used to select indicators.
}
\usage{
coverage(x, selection=NULL) 
plotcoverage(x, by=0.05, type="stat", max.order=NULL, add=FALSE , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class '\code{\link{indicators}}'.}
  \item{selection}{A logical vector restricting the set of indicators used to calculate the coverage.}
  \item{by}{Rate of increase in the predictive value threshold (At).}
  \item{type}{Specifies how to calculated the coverage: from the A statistic (\code{type = "stat"}), the lower bound of its confidence interval (\code{type = "lowerCI"}) or the upper bound (\code{type = "upperCI"}).}
  \item{max.order}{The maximum number of species conforming species combinations (allows examining the effects of increasing the order of combinations).}
  \item{add}{Flag indicating whether the function should draw on the existing plot.}
  \item{...}{Additional plotting parameters that are passed to the \code{plot} function.}
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
\seealso{\code{\link{indicators}}, \code{\link{pruneindicators}}}
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
sc= indicators(X=wetland[,sel], cluster=wetkm$cluster, group=1, verbose=TRUE, At=0.5, Bt=0.2)

## Determine the coverage of the selected set of indicators
coverage(sc)

## Plot the coverage against the threshold At
plotcoverage(sc)
plotcoverage(sc, max.order=2, add=TRUE, lty=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.