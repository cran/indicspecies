## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(indicspecies)
data(pigeons)

## -----------------------------------------------------------------------------
ls()

## -----------------------------------------------------------------------------
dfood

## ----fig=TRUE, fig.height=5, fig.width=5--------------------------------------
plot(hclust(dfood, method="average"), h=-1, xlab="", 
     ylab="Distance", main="", sub="", ylim=c(0.2,1))


## -----------------------------------------------------------------------------
diet.pop.barcelona <- colSums(diet.barcelona)
round(diet.pop.barcelona/sum(diet.pop.barcelona), dig=3)
diet.pop.moia <- colSums(diet.moia)
round(diet.pop.moia/sum(diet.pop.moia), dig=3)

## -----------------------------------------------------------------------------
nichevar(P=diet.barcelona, mode="single")
nichevar(P=diet.moia, mode="single")

## -----------------------------------------------------------------------------
popvar.barcelona <- nichevar(P=diet.barcelona, D=dfood, 
                            mode="single")
popvar.barcelona
popvar.moia <- nichevar(P=diet.moia, D=dfood, mode="single")
popvar.moia

## -----------------------------------------------------------------------------
nicheoverlap(P1=diet.barcelona, P2=diet.moia, mode="single")
nicheoverlap(P1=diet.barcelona, P2=diet.moia, mode="single", D = dfood)

## -----------------------------------------------------------------------------
round(sweep(diet.barcelona, 1, FUN="/", 
            rowSums(diet.barcelona)), dig=3)

## -----------------------------------------------------------------------------
round(sweep(diet.moia, 1, FUN="/", 
            rowSums(diet.moia)), dig=3)

## -----------------------------------------------------------------------------
indvar.barcelona <- nichevar(P=diet.barcelona, D=dfood)
summary(indvar.barcelona)
indvar.moia <- nichevar(P=diet.moia, D=dfood)
summary(indvar.moia)

## -----------------------------------------------------------------------------
wilcox.test(indvar.barcelona$B, indvar.moia$B)

## -----------------------------------------------------------------------------
Spec.barcelona <- mean(indvar.barcelona$B)/popvar.barcelona$B
Spec.barcelona
Spec.moia <- mean(indvar.moia$B)/popvar.moia$B
Spec.moia

## -----------------------------------------------------------------------------
Spec.ind.barcelona <- indvar.barcelona$B/popvar.barcelona$B
Spec.ind.moia <- indvar.moia$B/popvar.moia$B

## -----------------------------------------------------------------------------
wilcox.test(Spec.ind.barcelona, Spec.ind.moia)

## -----------------------------------------------------------------------------
O.barcelona <- nicheoverlap(diet.barcelona,D=dfood, mode="pairwise")
O.moia <- nicheoverlap(diet.moia,D=dfood, mode="pairwise")

## -----------------------------------------------------------------------------
mean(O.barcelona[lower.tri(O.barcelona)])
mean(O.moia[lower.tri(O.moia)])

## -----------------------------------------------------------------------------
O.barcelona.ind <- (rowSums(O.barcelona)-1)/(nrow(O.barcelona)-1)
summary(O.barcelona.ind)
O.moia.ind <- (rowSums(O.moia)-1)/(nrow(O.moia)-1)
summary(O.moia.ind)

## -----------------------------------------------------------------------------
wilcox.test(O.barcelona.ind, O.moia.ind)

