## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE,
    dev="jpeg", dpi = 72, fig.width = 4.5, fig.height = 3.5)
library(BiocStyle)

## ---- eval=FALSE--------------------------------------------------------------
#  # Install BiocManager if needed
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  # Install dittoSeq
#  BiocManager::install("dittoSeq")

## -----------------------------------------------------------------------------
## Download Data
library(scRNAseq)
sce <- BaronPancreasData()
# Trim to only 5 of the cell types for simplicity of vignette
sce <- sce[,sce$label %in% c(
    "acinar", "beta", "gamma", "delta", "ductal")]

## -----------------------------------------------------------------------------
## Some Quick Pre-processing
# Normalization.
library(scater)
sce <- logNormCounts(sce)

# Feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# PCA & UMAP
library(scater)
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)
sce <- runUMAP(sce, pca = 10)

# Clustering.
library(bluster)
sce$cluster <- clusterCells(sce, use.dimred='PCA',
    BLUSPARAM=NNGraphParam(cluster.fun="louvain"))

# Add some metadata common to Seurat objects
sce$nCount_RNA <- colSums(counts(sce))
sce$nFeature_RNA <- colSums(counts(sce)>0)
sce$percent.mito <- colSums(counts(sce)[grep("^MT-", rownames(sce)),])/sce$nCount_RNA 

sce

## -----------------------------------------------------------------------------
library(dittoSeq)
dittoDimPlot(sce, "donor")
dittoPlot(sce, "ENO1", group.by = "label")
dittoBarPlot(sce, "label", group.by = "donor")

## -----------------------------------------------------------------------------
# First, we'll just make some mock expression and conditions data
exp <- matrix(rpois(20000, 5), ncol=20)
colnames(exp) <- paste0("donor", seq_len(ncol(exp)))
rownames(exp) <- paste0("gene", seq_len(nrow(exp)))
logexp <- logexp <- log2(exp + 1)

pca <- matrix(rnorm(20000), nrow=20)

conditions <- factor(rep(1:4, 5))
sex <- c(rep("M", 9), rep("F", 11))

## -----------------------------------------------------------------------------
library(SummarizedExperiment)
bulkSE <- SummarizedExperiment(
    assays = list(counts = exp,
             logcounts = logexp),
    colData = data.frame(conditions = conditions,
                          sex = sex)
)

## -----------------------------------------------------------------------------
# dittoSeq import which allows
bulkSCE <- importDittoBulk(
    # x can be a DGEList, a DESeqDataSet, a SummarizedExperiment,
    #   or a list of data matrices
    x = list(counts = exp,
             logcounts = logexp),
    # Optional inputs:
    #   For adding metadata
    metadata = data.frame(conditions = conditions,
                          sex = sex),
    #   For adding dimensionality reductions
    reductions = list(pca = pca)
    )

## -----------------------------------------------------------------------------
# Add metadata (metadata can alternatively be added in this way)
bulkSCE$conditions <- conditions
bulkSCE$sex <- sex

# Add dimensionality reductions (can alternatively be added this way)
bulkSCE <- addDimReduction(
    object = bulkSCE,
    # (We aren't actually calculating PCA here.)
    embeddings = pca,
    name = "pca",
    key = "PC")

## -----------------------------------------------------------------------------
library(dittoSeq)
dittoDimPlot(bulkSCE, "sex", size = 3, do.ellipse = TRUE)
dittoBarPlot(bulkSCE, "sex", group.by = "conditions")
dittoBoxPlot(bulkSCE, "gene1", group.by = "sex")
dittoHeatmap(bulkSCE, getGenes(bulkSCE)[1:10],
    annot.by = c("conditions", "sex"))

## ---- eval = FALSE------------------------------------------------------------
#  # SummarizedExperiment dim-plots:
#  dittoDimPlot(
#      bulkSE,"sex", size = 3, do.ellipse = TRUE,
#      reduction.use = pca
#      )

## -----------------------------------------------------------------------------
# Retrieve all metadata slot names
getMetas(sce)
# Query for the presence of a metadata slot
isMeta("nCount_RNA", sce)
# Retrieve metadata values:
meta("label", sce)[1:10]
# Retrieve unique values of a metadata
metaLevels("label", sce)

## -----------------------------------------------------------------------------
# Retrieve all gene names
getGenes(sce)[1:10]
# Query for the presence of a gene(s)
isGene("CD3E", sce)
isGene(c("CD3E","ENO1","INS","non-gene"), sce, return.values = TRUE)
# Retrieve gene expression values:
gene("ENO1", sce)[1:10]

## -----------------------------------------------------------------------------
# Retrieve all dimensionality reductions
getReductions(sce)

## -----------------------------------------------------------------------------
# Getter
isBulk(sce)
isBulk(bulkSCE)

# Setter
mock_bulk <- setBulk(sce) # to bulk
isBulk(sce)
mock_sc <- setBulk(bulkSCE, set = FALSE) # to single-cell
isBulk(bulkSCE)

## ---- results = "hold"--------------------------------------------------------
dittoDimPlot(sce, "label", reduction.use = "PCA")
dittoDimPlot(sce, "ENO1")

## ---- results = "hold"--------------------------------------------------------
dittoScatterPlot(
    object = sce,
    x.var = "PPY", y.var = "INS",
    color.var = "label")
dittoScatterPlot(
    object = sce,
    x.var = "nCount_RNA", y.var = "nFeature_RNA",
    color.var = "percent.mito")

## -----------------------------------------------------------------------------
dittoDimPlot(sce, "cluster",
             
             do.label = TRUE,
             labels.repel = FALSE,
             
             add.trajectory.lineages = list(
                 c("9","3"),
                 c("8","7","2","4"),
                 c("8","7","1"),
                 c("5","11","6"),
                 c("10","0")),
             trajectory.cluster.meta = "cluster")

## ---- results = "hold"--------------------------------------------------------
dittoDimHex(sce)
dittoScatterHex(sce,
    x.var = "PPY", y.var = "INS")

## ---- results = "hold"--------------------------------------------------------
dittoDimHex(sce, "INS")
dittoScatterHex(
    object = sce,
    x.var = "PPY", y.var = "INS",
    color.var = "label",
    colors = c(1:4,7), max.density = 15)

## ---- results = "hold"--------------------------------------------------------
dittoPlot(sce, "ENO1", group.by = "label",
    plots = c("vlnplot", "jitter"))
dittoRidgePlot(sce, "ENO1", group.by = "label")
dittoBoxPlot(sce, "ENO1", group.by = "label")

## -----------------------------------------------------------------------------
dittoPlot(sce, "ENO1", group.by = "label",
    plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
    
    # change the color and size of jitter points
    jitter.color = "blue", jitter.size = 0.5,
    
    # change the outline color and width, and remove the fill of boxplots
    boxplot.color = "white", boxplot.width = 0.1,
    boxplot.fill = FALSE,
    
    # change how the violin plot widths are normalized across groups
    vlnplot.scaling = "count"
    )

## ---- results = "hold"--------------------------------------------------------
# dittoBarPlot
dittoBarPlot(sce, "label", group.by = "donor")
dittoBarPlot(sce, "label", group.by = "donor",
    scale = "count")

## ---- results = "hold"--------------------------------------------------------
# dittoFreqPlot
sce$mock.donor.group <- ifelse(sce$donor %in% unique(sce$donor)[1:2], "A", "B")
dittoFreqPlot(sce, "label",
    sample.by = "donor", group.by = "mock.donor.group")

## ---- results = "hold"--------------------------------------------------------
# Pick Genes
genes <- c("SST", "REG1A", "PPY", "INS", "CELA3A", "PRSS2", "CTRB1",
    "CPA1", "CTRB2" , "REG3A", "REG1B", "PRSS1", "GCG", "CPB1",
    "SPINK1", "CELA3B", "CLPS", "OLFM4", "ACTG1", "FTL")

# Annotating and ordering cells by some meaningful feature(s):
dittoHeatmap(sce, genes,
    annot.by = c("label", "donor"))
dittoHeatmap(sce, genes,
    annot.by = c("label", "donor"),
    order.by = "donor")

## -----------------------------------------------------------------------------
# Add annotations
dittoHeatmap(sce, genes,
    annot.by = c("label", "donor"),
    scaled.to.max = TRUE,
    show_colnames = FALSE,
    show_rownames = FALSE)

## -----------------------------------------------------------------------------
# Highlight certain genes
dittoHeatmap(sce, genes, annot.by = c("label", "donor"),
    highlight.features = genes[1:3],
    complex = TRUE)

## -----------------------------------------------------------------------------
# seurat <- as.Seurat(sce)
# Idents(seurat) <- "label"
# delta.marker.table <- FindMarkers(seurat, ident.1 = "delta")
# delta.genes <- rownames(delta.marker.table)[1:20]
# Idents(seurat) <- "seurat_clusters"

delta.genes <- c(
    "SST", "RBP4", "LEPR", "PAPPA2", "LY6H",
    "CBLN4", "GPX3", "BCHE", "HHEX", "DPYSL3",
    "SERPINA1", "SEC11C", "ANXA2", "CHGB", "RGS2",
    "FXYD6", "KCNIP1", "SMOC1", "RPL10", "LRFN5")

## -----------------------------------------------------------------------------
dittoDotPlot(sce, vars = delta.genes, group.by = "label")
dittoDotPlot(sce, vars = delta.genes, group.by = "label",
    scale = FALSE)

## -----------------------------------------------------------------------------
multi_dittoPlot(sce, delta.genes[1:6], group.by = "label",
    vlnplot.lineweight = 0.2, jitter.size = 0.3)
dittoPlotVarsAcrossGroups(sce, delta.genes, group.by = "label",
    main = "Delta-cell Markers")

## ---- results = "hold"--------------------------------------------------------
multi_dittoDimPlot(sce, delta.genes[1:6])
multi_dittoDimPlotVaryCells(sce, delta.genes[1],
    vary.cells.meta = "label")

## -----------------------------------------------------------------------------
# Original
dittoBarPlot(sce, "label", group.by = "donor", scale = "count")

# First 10 cells
dittoBarPlot(sce, "label", group.by = "donor", scale = "count",
    # String method
    cells.use = colnames(sce)[1:10]
    # Index method, which would achieve the same effect
    # cells.use = 1:10
    )

# Acinar cells only
dittoBarPlot(sce, "label", group.by = "donor", scale = "count",
    # Logical method
    cells.use = meta("label", sce) == "acinar")

## -----------------------------------------------------------------------------
dittoPlot(sce, "PPY", group.by = "donor", 
    split.by = "label")
dittoDimPlot(sce, "PPY",
    split.by = c("donor", "label"))

## -----------------------------------------------------------------------------
dittoPlot(sce, "PPY", group.by = "donor", 
    split.by = "label",
    split.adjust = list(scales = "free_y"), max = NA)

## ---- fig.height=7------------------------------------------------------------
dittoRidgePlot(sce, "PPY", group.by = "donor", 
    split.by = "label",
    split.ncol = 1)

## -----------------------------------------------------------------------------
dittoBarPlot(sce, "label", group.by = "donor",
    main = "Encounters",
    sub = "By Type",
    xlab = NULL, # NULL = remove
    ylab = "Generation 1",
    x.labels = c("Ash", "Misty", "Jessie", "James"),
    legend.title = "Types",
    var.labels.rename = c("Fire", "Water", "Grass", "Electric", "Psychic"),
    x.labels.rotate = FALSE)

## ---- results="hold"----------------------------------------------------------
# original - discrete
dittoDimPlot(sce, "label")
# swapped colors
dittoDimPlot(sce, "label",
    colors = 5:1)
# different colors
dittoDimPlot(sce, "label",
    color.panel = c("red", "orange", "purple", "yellow", "skyblue"))

## ---- results="hold"----------------------------------------------------------
# original - expression
dittoDimPlot(sce, "INS")
# different colors
dittoDimPlot(sce, "INS",
    max.color = "red", min.color = "gray90")

## -----------------------------------------------------------------------------
dittoBarPlot(sce, "label", group.by = "donor",
    data.out = TRUE)

## -----------------------------------------------------------------------------
dittoHeatmap(sce, c("SST","CPE","GPX3"), cells.use = colnames(sce)[1:5],
    data.out = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  # These can be finicky to render in knitting, but still, example code:
#  dittoDimPlot(sce, "INS",
#      do.hover = TRUE,
#      hover.data = c("label", "donor", "ENO1", "cluster", "nCount_RNA"))
#  dittoBarPlot(sce, "label", group.by = "donor",
#      do.hover = TRUE)

## -----------------------------------------------------------------------------
# Note: dpi gets re-set by the styling code of this vignette, so this is
#   just a code example, but the plot won't be quite matched.
dittoDimPlot(sce, "label",
    do.raster = TRUE,
    raster.dpi = 300)

## -----------------------------------------------------------------------------
dittoHeatmap(sce, genes, scaled.to.max = TRUE,
    complex = TRUE,
    use_raster = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()

