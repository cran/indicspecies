\encoding{UTF-8}
\name{summary.multipatt}
\alias{summary.multipatt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing multi-level pattern analysis}
\description{
 This is a \code{summary} function for \code{multipatt} objects. 
}
\usage{
\method{summary}{multipatt}(object, alpha=0.05, minstat=0.0, indvalcomp=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned by function \code{\link{multipatt}}.}
  \item{alpha}{Significance level for selecting species in the summary.}
  \item{minstat}{Minimum value of the statistic for selecting species in the summary.}
  \item{indvalcomp}{Flag used to report components "A" (or "A.g") and "B" of indicator values. This argument is effective only if "multipatt" function was executed using \code{func = "IndVal"} or \code{func = "IndVal.g"}.}
  \item{...}{Additional arguments of the \code{summary} function.}
}
\details{
   The summary function for \code{multipatt} objects presents the results of the analysis by grouping species following their associated site group combination. Only significantly associated species are listed for each combination. The user can limit the species shown by setting a significance level and the minimum value for the strength of the association. 
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Biodiversity and Landscape Ecology Laboratory, Forest Science Center of Catalonia
}
\seealso{\code{\link{strassoc}}, \code{\link{signassoc}}, \code{\link{multipatt}}}
\examples{
library(stats)

data(wetland) ## Loads species data

wetkm = kmeans(wetland, centers=3) ## Creates three clusters using kmeans

wetpt = multipatt(wetland, wetkm$cluster) ## Runs the combination analysis using IndVal.g as statistic

summary(wetpt) ## Lists those species with significant association to one combination

summary(wetpt, indvalcomp=TRUE) ## Lists those species with significant association to one combination. Includes indval components
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.