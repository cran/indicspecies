## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(indicspecies)

## -----------------------------------------------------------------------------
data(wetland)

## -----------------------------------------------------------------------------
groups <- c(rep(1, 17), rep(2, 14), rep(3,10))
groups

## -----------------------------------------------------------------------------
wetkm <- kmeans(wetland, centers=3)
groupskm <- wetkm$cluster
groupskm

## -----------------------------------------------------------------------------
indval <- multipatt(wetland, groups, 
                    control = how(nperm=999)) 

## -----------------------------------------------------------------------------
summary(indval) 

## -----------------------------------------------------------------------------
summary(indval, indvalcomp=TRUE)

## -----------------------------------------------------------------------------
summary(indval, alpha=1)

## -----------------------------------------------------------------------------
indval$sign

## -----------------------------------------------------------------------------
wetlandpa <- ifelse(wetland>0,1,0)
phi <- multipatt(wetlandpa, groups, func = "r", 
                 control = how(nperm=999)) 

## -----------------------------------------------------------------------------
phi <- multipatt(wetlandpa, groups, func = "r.g", 
                 control = how(nperm=999)) 

## -----------------------------------------------------------------------------
summary(phi)

## -----------------------------------------------------------------------------
round(head(phi$str),3)

## -----------------------------------------------------------------------------
round(head(indval$str),3)

## -----------------------------------------------------------------------------
phi <- multipatt(wetlandpa, groups, func = "r.g", 
                 control = how(nperm=999), 
                 allow.negative = TRUE) 

## -----------------------------------------------------------------------------
summary(phi)

## -----------------------------------------------------------------------------
indvalori <- multipatt(wetland, groups, duleg = TRUE, 
                       control = how(nperm=999)) 
summary(indvalori)

## -----------------------------------------------------------------------------
indvalrest <- multipatt(wetland, groups, max.order = 2, 
                        control = how(nperm=999)) 
summary(indvalrest)

## -----------------------------------------------------------------------------
indvalrest <- multipatt(wetland, groups, restcomb = c(1,2,3,5,6), 
                        control = how(nperm=999)) 
summary(indvalrest)

## -----------------------------------------------------------------------------
indvalrest$sign

## -----------------------------------------------------------------------------
prefstat <- strassoc(wetland, cluster=groups, func="A.g")
round(head(prefstat),3)

## -----------------------------------------------------------------------------
prefstat <- strassoc(wetland, cluster=groups, func="A.g", 
                     nboot.ci = 199)
round(head(prefstat$lowerCI),3)
round(head(prefstat$upperCI),3)

## -----------------------------------------------------------------------------
prefsign <- signassoc(wetland, cluster=groups,  alternative = "two.sided", 
                      control = how(nperm=199)) 
head(prefsign)

## -----------------------------------------------------------------------------
coverage(wetland, indvalori)

## -----------------------------------------------------------------------------
coverage(wetland, indvalori, At = 0.8, alpha = 0.05)

## ----fig = TRUE, fig.width = 5, fig.height = 5--------------------------------
plotcoverage(x=wetland, y=indvalori, group="1", lty=1)
plotcoverage(x=wetland, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=wetland, y=indvalori, group="3", lty=3, col="red", add=TRUE)
legend(x = 0.1, y=30, 
       legend=c("group 1","group 2", "group 3"),
       lty=c(1,2,3), col=c("black","blue","red"), bty="n")

## -----------------------------------------------------------------------------
wetcomb <- combinespecies(wetland, max.order = 2)$XC
dim(wetcomb)

## -----------------------------------------------------------------------------
indvalspcomb <- multipatt(wetcomb, groups, duleg = TRUE, 
                          control = how(nperm=999))
summary(indvalspcomb, indvalcomp = TRUE)

## -----------------------------------------------------------------------------
sc <- indicators(X=wetland, cluster=groups, group=2, 
                 max.order = 3, verbose=TRUE, 
                 At=0.5, Bt=0.2)

## -----------------------------------------------------------------------------
print(sc, sqrtIVt = 0.6)

## -----------------------------------------------------------------------------
coverage(sc)

## -----------------------------------------------------------------------------
coverage(sc, At=0.8, alpha =0.05)

## ----fig = TRUE, fig.width = 5, fig.height = 5--------------------------------
plotcoverage(sc)
plotcoverage(sc, max.order=1, add=TRUE, lty=2, col="red")
legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"), 
       lty=c(1,2), col=c("black","red"), bty="n")

## -----------------------------------------------------------------------------
sc2 <- pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)

## -----------------------------------------------------------------------------
p <- predict(sc2, wetland)

## -----------------------------------------------------------------------------
p <- predict(sc2)

## -----------------------------------------------------------------------------
pcv <- predict(sc2, cv=TRUE)

## -----------------------------------------------------------------------------
data.frame(Group2 = as.numeric(wetkm$cluster==2), Prob = p, Prob_CV = pcv)

