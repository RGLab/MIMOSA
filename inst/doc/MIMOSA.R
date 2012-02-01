### R code from vignette source 'MIMOSA.Rnw'

###################################################
### code chunk number 1: loaddata
###################################################
require(MIMOSA)
head(HVTN049)


###################################################
### code chunk number 2: reshape
###################################################
HVTN049.ics<-ICS(HVTN049)
show(HVTN049.ics)


###################################################
### code chunk number 3: extract
###################################################
data<-flowContrasts(HVTN049.ics, control="negctrl1",stim="ENV-1-PTEG",subset=c("cd4","IL2"))
head(data)


###################################################
### code chunk number 4: fit
###################################################
result<-BetaMix(data,scl=100)
show(result)


###################################################
### code chunk number 5: fdr
###################################################
table(result@fdr<0.01)


###################################################
### code chunk number 6: volcano
###################################################
res<-volcanoPlot(result)


###################################################
### code chunk number 7: densityplots
###################################################
densityplot(result,c(24,1))


###################################################
### code chunk number 8: contours
###################################################
contourPlot(result,which=24,nsamps=2000,n=30)


###################################################
### code chunk number 9: polyfunctionality
###################################################
poly<-testPolyfunctionality(ics=HVTN049.ics,stim="ENV-1-PTEG",control="negctrl1",cytokineA="IL2",cytokineB="IFNg",subset=c("cd4"),or="IFNg.IL2")


###################################################
### code chunk number 10: unifunctionality
###################################################
ifng<-BetaMix(flowContrasts(HVTN049.ics,stim="ENV-1-PTEG",control="negctrl1",subset=c("cd4","IFNg")),scl=100)
il2<-BetaMix(flowContrasts(HVTN049.ics,stim="ENV-1-PTEG",control="negctrl1",subset=c("cd4","IL2")),scl=100)


###################################################
### code chunk number 11: summary
###################################################
summary(poly)
summary(ifng)
summary(il2)


