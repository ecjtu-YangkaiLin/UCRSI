library(Seurat)
library(ggpubr)
library(aricode)
library(Matrix)
library(reticulate)
library(umap)

setwd("~/Project/UCRSI")
source(paste0(getwd(), "/FeatureSelection.R"))
source(paste0(getwd(), "/UCRSI.R"))

# Parameters -------------------------------------------------------------------

# Change this to your own python environment path, requires sklearn and pyamg packages
python_env = "/home/linyangkai/virtualenv/venv-single-cell" 
n_cluster = 9 # Number of clusters
n_neighbors = 10
min_Gene = 4
output_distance = T # Output two distance lists that can be used to improve visualization and measure tumor heterogeneity
alpha_start = 0
alpha_end = 1
alpha_step = 0.1
nJobs = 31 # How many cores are used for distance calculation



# Download dataset (GSE103322)--------------------------------------------------
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103322&format=file&file=GSE103322%5FHNSCC%5Fall%5Fdata%2Etxt%2Egz"
data.file = paste0(getwd(), "/GSE103322.gz")
download.file(url, data.file)
counts = read.table(gzfile(data.file), sep = "\t", header = T, row.names = 1)


# Creating a seurat object -----------------------------------------------------
counts = counts[, counts[3, ] != 0 | counts[4, ] != 0]
celltype = as.vector(as.matrix(counts[5, ]))
celltype[celltype == "0"] = "Cancer"
celltype[celltype == "-Fibroblast"] = "Fibroblast"
counts = counts[-1:-5, ]
seurat.obj = CreateSeuratObject(
  counts = counts, 
  project = "GSE103322: Single cell RNA-seq analysis of head and neck cancer")
seurat.obj@meta.data[["orig.ident"]] = as.factor(celltype)
seurat.obj = NormalizeData(seurat.obj)


# Distance calculation ---------------------------------------------------------
features = features_selection(seurat.obj, nbin = max(10, n_cluster * 2))
UCRSI.results = UCRSI(as.matrix(seurat.obj@assays[["RNA"]]@data[features$feature.genes, ]),
                      minGene = min_Gene,
                      n_neighbors = n_neighbors,
                      outputDistance = output_distance, 
                      alpha_start = alpha_start, 
                      alpha_end = alpha_end, 
                      alpha_step = alpha_step, 
                      nJobs = nJobs,
                      infValue = 0x7fffffff)


# Clustering using sklearn -----------------------------------------------------
if(!py_available()) {
  use_virtualenv(python_env, required = T)
  py_available(T)
}
pyamg = import("pyamg")
sklearn.cluster = import("sklearn.cluster")
seurat.obj@meta.data[["UCRSI"]] = as.factor(sklearn.cluster$SpectralClustering(
  n_clusters=as.integer(n_cluster), 
  n_neighbors = as.integer(n_neighbors),
  affinity="precomputed",
  assign_labels='discretize', 
  n_jobs = as.integer(nJobs),
  eigen_solver = "amg",
  random_state=as.integer(1))$fit_predict(UCRSI.results$adjacency))


# Performance evaluation -------------------------------------------------------
NMI(seurat.obj@meta.data[["orig.ident"]], seurat.obj@meta.data[["UCRSI"]], variant = "sum")
AMI(seurat.obj@meta.data[["orig.ident"]], seurat.obj@meta.data[["UCRSI"]])
ARI(seurat.obj@meta.data[["orig.ident"]], seurat.obj@meta.data[["UCRSI"]])



# Visualization using features genes -------------------------------------------
seurat.obj = RunUMAP(seurat.obj, features = features$feature.genes)
ggarrange(DimPlot(seurat.obj, reduction = "umap", label = T, pt.size = 0.5, group.by="orig.ident"),
          DimPlot(seurat.obj, reduction = "umap", label = T, pt.size = 0.5, group.by="UCRSI"),
          ncol = 2, nrow = 1)


# Visualization using UCRSI distance -------------------------------------------
custom.config = umap.defaults
custom.config$n_neighbors = 30
custom.config$min_dist = 0.3

UCRSI.1 = distanceListToMatrix(UCRSI.results[["distance1"]]) # Convert a list of distances into a complete matrix
UCRSI.2 = distanceListToMatrix(UCRSI.results[["distance2"]])

USRCI.umap.1 = as.data.frame(umap(UCRSI.1, input="dist", config=custom.config)[["layout"]])
seurat.obj@reductions[["USRCI.umap.1"]] = seurat.obj@reductions[["umap"]]
seurat.obj@reductions[["USRCI.umap.1"]]@cell.embeddings[, 1] = USRCI.umap.1$V1
seurat.obj@reductions[["USRCI.umap.1"]]@cell.embeddings[, 2] = USRCI.umap.1$V2
p1 = DimPlot(seurat.obj, reduction = "USRCI.umap.1", label = T, pt.size = 0.5, group.by="orig.ident") + ggtitle(bquote('UMAP'~(UCRSI[1])))


USRCI.umap.2 = as.data.frame(umap(UCRSI.2, input="dist", config=custom.config)[["layout"]])
seurat.obj@reductions[["USRCI.umap.2"]] = seurat.obj@reductions[["umap"]]
seurat.obj@reductions[["USRCI.umap.2"]]@cell.embeddings[, 1] = USRCI.umap.2$V1
seurat.obj@reductions[["USRCI.umap.2"]]@cell.embeddings[, 2] = USRCI.umap.2$V2
p2 = DimPlot(seurat.obj, reduction = "USRCI.umap.2", label = T, pt.size = 0.5, group.by="orig.ident") + ggtitle(bquote('UMAP'~(UCRSI[2])))


USRCI.umap = as.data.frame(umap(UCRSI.1+UCRSI.2, input="dist", config=custom.config)[["layout"]])
seurat.obj@reductions[["USRCI.umap"]] = seurat.obj@reductions[["umap"]]
seurat.obj@reductions[["USRCI.umap"]]@cell.embeddings[, 1] = USRCI.umap$V1
seurat.obj@reductions[["USRCI.umap"]]@cell.embeddings[, 2] = USRCI.umap$V2
p3 = DimPlot(seurat.obj, reduction = "USRCI.umap", label = T, pt.size = 0.5, group.by="orig.ident") + ggtitle(bquote('UMAP'~(UCRSI[1]+UCRSI[2])))

ggarrange(p1, p2, p3, labels = "AUTO")
