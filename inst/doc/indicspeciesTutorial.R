### R code from vignette source 'indicspeciesTutorial.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: indicspeciesTutorial.Rnw:18-19
###################################################
options(width=67)


###################################################
### code chunk number 2: indicspeciesTutorial.Rnw:27-28
###################################################
library(indicspecies)


###################################################
### code chunk number 3: indicspeciesTutorial.Rnw:35-36
###################################################
data(wetland)


###################################################
### code chunk number 4: indicspeciesTutorial.Rnw:42-44
###################################################
groups = c(rep(1, 17), rep(2, 14), rep(3,10))
groups


###################################################
### code chunk number 5: indicspeciesTutorial.Rnw:47-50
###################################################
wetkm = kmeans(wetland, centers=3)
groupskm = wetkm$cluster
groupskm


###################################################
### code chunk number 6: indicspeciesTutorial.Rnw:61-63
###################################################
indval = multipatt(wetland, groups, 
                   control = how(nperm=999)) 


###################################################
### code chunk number 7: indicspeciesTutorial.Rnw:70-71
###################################################
summary(indval) 


###################################################
### code chunk number 8: indicspeciesTutorial.Rnw:79-80
###################################################
summary(indval, indvalcomp=TRUE)


###################################################
### code chunk number 9: indicspeciesTutorial.Rnw:86-87
###################################################
summary(indval, alpha=1)


###################################################
### code chunk number 10: indicspeciesTutorial.Rnw:90-91
###################################################
indval$sign


###################################################
### code chunk number 11: indicspeciesTutorial.Rnw:97-100
###################################################
wetlandpa = ifelse(wetland>0,1,0)
phi = multipatt(wetlandpa, groups, func = "r", 
                control = how(nperm=999)) 


###################################################
### code chunk number 12: indicspeciesTutorial.Rnw:105-107
###################################################
phi = multipatt(wetlandpa, groups, func = "r.g", 
                control = how(nperm=999)) 


###################################################
### code chunk number 13: indicspeciesTutorial.Rnw:113-114
###################################################
summary(phi)


###################################################
### code chunk number 14: indicspeciesTutorial.Rnw:119-120
###################################################
round(head(phi$str),3)


###################################################
### code chunk number 15: indicspeciesTutorial.Rnw:123-124
###################################################
round(head(indval$str),3)


###################################################
### code chunk number 16: indicspeciesTutorial.Rnw:134-137
###################################################
indvalori = multipatt(wetland, groups, duleg = TRUE, 
                      control = how(nperm=999)) 
summary(indvalori)


###################################################
### code chunk number 17: indicspeciesTutorial.Rnw:142-145
###################################################
indvalrest = multipatt(wetland, groups, max.order = 2, 
                       control = how(nperm=999)) 
summary(indvalrest)


###################################################
### code chunk number 18: indicspeciesTutorial.Rnw:153-156
###################################################
indvalrest = multipatt(wetland, groups, restcomb = c(1,2,3,5,6), 
                       control = how(nperm=999)) 
summary(indvalrest)


###################################################
### code chunk number 19: indicspeciesTutorial.Rnw:159-160
###################################################
indvalrest$sign


###################################################
### code chunk number 20: indicspeciesTutorial.Rnw:168-170
###################################################
prefstat = strassoc(wetland, cluster=groups, func="A.g")
round(head(prefstat),3)


###################################################
### code chunk number 21: indicspeciesTutorial.Rnw:173-177
###################################################
prefstat = strassoc(wetland, cluster=groups, func="A.g", 
                    nboot.ci = 199)
round(head(prefstat$lowerCI),3)
round(head(prefstat$upperCI),3)


###################################################
### code chunk number 22: indicspeciesTutorial.Rnw:184-187
###################################################
prefsign = signassoc(wetland, cluster=groups,  alternative = "two.sided", 
                     control = how(nperm=199)) 
head(prefsign)


###################################################
### code chunk number 23: indicspeciesTutorial.Rnw:195-196
###################################################
coverage(wetland, indvalori)


###################################################
### code chunk number 24: indicspeciesTutorial.Rnw:201-202
###################################################
coverage(wetland, indvalori, At = 0.8, alpha = 0.05)


###################################################
### code chunk number 25: indicspeciesTutorial.Rnw:210-216
###################################################
plotcoverage(x=wetland, y=indvalori, group="1", lty=1)
plotcoverage(x=wetland, y=indvalori, group="2", lty=2, col="blue", add=TRUE)
plotcoverage(x=wetland, y=indvalori, group="3", lty=3, col="red", add=TRUE)
legend(x = 0.01, y=20, 
       legend=c("group 1","group 2", "group 3"),
       lty=c(1,2,3), col=c("black","blue","red"), bty="n")


###################################################
### code chunk number 26: indicspeciesTutorial.Rnw:228-230
###################################################
wetcomb = combinespecies(wetland, max.order = 2)$XC
dim(wetcomb)


###################################################
### code chunk number 27: indicspeciesTutorial.Rnw:233-236
###################################################
indvalspcomb = multipatt(wetcomb, groups, duleg = TRUE, 
                       control = how(nperm=999))
summary(indvalspcomb, indvalcomp = TRUE)


###################################################
### code chunk number 28: indicspeciesTutorial.Rnw:242-245
###################################################
sc= indicators(X=wetland, cluster=groups, group=2, 
               max.order = 3, verbose=TRUE, 
               At=0.5, Bt=0.2)


###################################################
### code chunk number 29: indicspeciesTutorial.Rnw:248-249
###################################################
print(sc, sqrtIVt = 0.6)


###################################################
### code chunk number 30: indicspeciesTutorial.Rnw:257-258
###################################################
coverage(sc)


###################################################
### code chunk number 31: indicspeciesTutorial.Rnw:261-262
###################################################
coverage(sc, At=0.8, alpha =0.05)


###################################################
### code chunk number 32: indicspeciesTutorial.Rnw:268-272
###################################################
plotcoverage(sc)
plotcoverage(sc, max.order=1, add=TRUE, lty=2, col="red")
legend(x=0.1, y=20, legend=c("Species combinations","Species singletons"), 
       lty=c(1,2), col=c("black","red"), bty="n")


###################################################
### code chunk number 33: indicspeciesTutorial.Rnw:279-281
###################################################
sc2=pruneindicators(sc, At=0.8, Bt=0.2, verbose=TRUE)
print(sc2)


###################################################
### code chunk number 34: indicspeciesTutorial.Rnw:287-288
###################################################
p<-predict(sc2, wetland)


###################################################
### code chunk number 35: indicspeciesTutorial.Rnw:291-292
###################################################
p<-predict(sc2)


###################################################
### code chunk number 36: indicspeciesTutorial.Rnw:295-296
###################################################
pcv<-predict(sc2, cv=TRUE)


###################################################
### code chunk number 37: indicspeciesTutorial.Rnw:299-300
###################################################
data.frame(Group2 = as.numeric(wetkm$cluster==2), Prob = p, Prob_CV = pcv)


