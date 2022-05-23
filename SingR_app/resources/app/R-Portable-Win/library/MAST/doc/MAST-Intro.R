## ----long-example, warning=FALSE, echo=FALSE----------------------------------
knitr::opts_chunk$set(error=FALSE, echo=TRUE)
suppressPackageStartupMessages({
library(MAST)
library(data.table)
})

## ----long-example2,warning=FALSE----------------------------------------------
data(vbeta)
colnames(vbeta)
vbeta <- computeEtFromCt(vbeta)
vbeta.fa <- FromFlatDF(vbeta, idvars=c("Subject.ID", "Chip.Number", "Well"),
                          primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                          geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                          phenovars=c('Stim.Condition','Time'), id='vbeta all', class='FluidigmAssay')
print(vbeta.fa)

## ----examineMeta--------------------------------------------------------------
head(rowData(vbeta.fa),3)
head(colData(vbeta.fa),3)

## ----subsets,warning=FALSE----------------------------------------------------
sub1 <- vbeta.fa[,1:10]
show(sub1)
sub2 <- subset(vbeta.fa, Well=='A01')
show(sub2)
sub3 <- vbeta.fa[6:10, 1:10]
show(sub3)

## ----split, warning=FALSE-----------------------------------------------------
sp1 <- split(vbeta.fa, 'Subject.ID')
show(sp1)

## ----combine,warning=FALSE,echo=FALSE-----------------------------------------
cbind(sp1[[1]],sp1[[2]])

## ----splitbyncells,warning=FALSE,fig.height=4, fig.width=4--------------------
vbeta.split<-split(vbeta.fa,"Number.of.Cells")
#see default parameters for plotSCAConcordance
plotSCAConcordance(vbeta.split[[1]],vbeta.split[[2]],
                   filterCriteria=list(nOutlier = 1, sigmaContinuous = 9,
                       sigmaProportion = 9))

## ----otherFiltering, warning=FALSE--------------------------------------------
vbeta.fa
## Split by 'ncells', apply to each component, then recombine
vbeta.filtered <- mast_filter(vbeta.fa, groups='ncells')
## Returned as boolean matrix
was.filtered <- mast_filter(vbeta.fa, apply_filter=FALSE)
## Wells filtered for being discrete outliers
head(subset(was.filtered, pctout))

## ----burdenOfFiltering,warning=FALSE,fig.width=4,fig.height=4-----------------
burdenOfFiltering(vbeta.fa, 'ncells', byGroup=TRUE)

## ----zlmArgs------------------------------------------------------------------
vbeta.1 <- subset(vbeta.fa, ncells==1)
## Consider the first 20 genes
vbeta.1 <- vbeta.1[1:20,] 
head(colData(vbeta.1))

## ----zlmExample, warning=FALSE, message=FALSE, fig.width=6, fig.height=6------
library(ggplot2)
zlm.output <- zlm(~ Population + Subject.ID, vbeta.1,)
show(zlm.output)
## returns a data.table with a summary of the fit
coefAndCI <- summary(zlm.output, logFC=FALSE)$datatable
coefAndCI <- coefAndCI[contrast != '(Intercept)',]
coefAndCI[,contrast:=abbreviate(contrast)]
ggplot(coefAndCI, aes(x=contrast, y=coef, ymin=ci.lo, ymax=ci.hi, col=component))+
    geom_pointrange(position=position_dodge(width=.5)) +facet_wrap(~primerid) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + coord_cartesian(ylim=c(-3, 3))

## ----tests, fig.width=4, fig.height=5-----------------------------------------
zlm.lr <- lrTest(zlm.output, 'Population')
dimnames(zlm.lr)
pvalue <- ggplot(melt(zlm.lr[,,'Pr(>Chisq)']), aes(x=primerid, y=-log10(value)))+
    geom_bar(stat='identity')+facet_wrap(~test.type) + coord_flip()
print(pvalue)

## ----lmerExample , warning=FALSE, message=FALSE, eval=FALSE-------------------
#  library(lme4)
#  lmer.output <- zlm(~ Stim.Condition +(1|Subject.ID), vbeta.1, method='glmer', ebayes=FALSE)

## ----advancedArmExample, fig.width=4, fig.height=4, message=FALSE, warning=FALSE----
 orig_results <- zlm(~Stim.Condition, vbeta.1)
 dp <- defaultPrior('Stim.ConditionUnstim')
 new_results <-  zlm(~Stim.Condition, vbeta.1, useContinuousBayes=TRUE,coefPrior=dp)
 qplot(x=coef(orig_results, 'C')[, 'Stim.ConditionUnstim'],
       y=coef(new_results, 'C')[, 'Stim.ConditionUnstim'],
       color=coef(new_results, 'D')[,'(Intercept)']) +
     xlab('Default prior') + ylab('Informative Prior') +
     geom_abline(slope=1, lty=2) + scale_color_continuous('Baseline\nlog odds\nof expression')

## ----LRTexample, eval=TRUE, error=TRUE----------------------------------------
two.sample <- LRT(vbeta.1, 'Population', referent='CD154+VbetaResponsive')
head(two.sample) 

