## ----init, results = 'hide'---------------------------------------------------
library(MAST)

knitr::opts_chunk$set(message = FALSE,error = FALSE,warning = FALSE)
data(maits, package='MAST')
unlog <- function(x) ceiling(2^x - 1)
sca_raw =  FromMatrix(t(maits$expressionmat), maits$cdat, maits$fdat)
assays(sca_raw)$counts = unlog(assay(sca_raw))
assayNames(sca_raw)

## ----scaterQC,results='hide'--------------------------------------------------
library(scater)
sca_raw = addPerCellQC(sca_raw)
plotColData(sca_raw, y="detected", x="total")

## -----------------------------------------------------------------------------
sca_raw = runPCA(sca_raw, ncomponents=5, exprs_values = 'et')
plotReducedDim(sca_raw, dimred = 'PCA', colour_by = 'condition')

## -----------------------------------------------------------------------------
  example_sce = mockSCE() 

example_sce = logNormCounts(example_sce)
sca = SceToSingleCellAssay(example_sce)

## -----------------------------------------------------------------------------
zlm( ~ Treatment, sca = sca, exprs_value = 'logcounts')

## -----------------------------------------------------------------------------
assayNames(sca)

## -----------------------------------------------------------------------------
library(Matrix)
sca_sparse = FromMatrix(
    exprsArray = list(et = Matrix(t(maits$expressionmat), sparse = TRUE)),
    cData = maits$cdat, fData = maits$fdat)
class(assay(sca_sparse))
regular_time = system.time(zlm( ~ condition, sca = sca_raw[1:100,]))
sparse_time = system.time(zlm( ~ condition, sca = sca_sparse[1:100,]))

## -----------------------------------------------------------------------------
library(DelayedArray)
library(HDF5Array)
hd5_dat = as(t(maits$expressionmat), 'HDF5Array')
DelayedArray::seed(hd5_dat)

## -----------------------------------------------------------------------------
sca_delay = FromMatrix(
    exprsArray = list(et = hd5_dat),
     cData = maits$cdat, fData = maits$fdat)
class(assay(sca_sparse))

hd5_time = system.time(zlm( ~ condition, sca = sca_delay[1:100,]))

## -----------------------------------------------------------------------------
knitr::kable(data.frame(method = c('Dense', 'Sparse', 'HDF5'), 'user time(s)' =c( regular_time[1], sparse_time[1], hd5_time[1]), check.names = FALSE))


## ----zinbwave-----------------------------------------------------------------
library(zinbwave)
feature_var = apply(assay(sca_raw), 1, var)
sca_top500 = sca_raw[rank(-feature_var)<=500,]
zw = zinbwave(Y = sca_top500, X = '~1', which_assay = 'counts', K = 2, normalizedValues = TRUE)


## ---- message=FALSE-----------------------------------------------------------
rd = data.frame(reducedDim(zw, 'PCA'), reducedDim(zw, 'zinbwave'), colData(zw))
GGally::ggpairs(rd, columns = c('PC1', 'PC2', 'W1', 'W2'), mapping = aes(color = condition))

## ---- results = 'hide'--------------------------------------------------------
colData(zw) = cbind(colData(zw), reducedDim(zw, 'zinbwave'))
zw = SceToSingleCellAssay(zw)
zz = zlm(~W1 + W2, sca = zw, exprs_values = 'et')

## ---- results = 'asis'--------------------------------------------------------
ss = summary(zz)
knitr::kable(print(ss))


## -----------------------------------------------------------------------------
library(dplyr)
library(data.table)
top5 = ss$datatable %>% filter(component=='logFC', contrast %like% 'W') %>% arrange(-abs(z)) %>% head(n=5) %>% left_join(rowData(zw) %>% as.data.table())
dat = zw[top5$primerid,] %>% as('data.table')
dat = dat[,!duplicated(colnames(dat)),with = FALSE]
plt = ggplot(dat, aes(x=W1, color = condition)) + geom_point() + facet_wrap(~symbolid)


## -----------------------------------------------------------------------------
plt + aes(y = et)

## -----------------------------------------------------------------------------
plt + aes(y = normalizedValues)

