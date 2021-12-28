features_selection = function(seurat_obj, min.mean = 0.1, nbin = 10) {
  seurat_obj = FindVariableFeatures(seurat_obj, mean.cutoff = c(0.1, 8), selection.method = "mvp", binning.method="equal_frequency")
  mvp.mean = seurat_obj@assays[["RNA"]]@meta.features[["mvp.mean"]]
  mvp.dispersion = seurat_obj@assays[["RNA"]]@meta.features[["mvp.dispersion"]]
  data = data.frame(x = mvp.mean, y = mvp.dispersion)
  figures = ggplot(data=data, aes(x=x, y=y)) + geom_point(size=0.5) + xlab("Average Expression") + ylab("Dispersion")
  
  gene_filter = which(mvp.mean > min.mean)
  
  mvp.mean = mvp.mean[gene_filter]
  mvp.dispersion = mvp.dispersion[gene_filter]
  
  gene.n = length(gene_filter)
  gene.list = rownames(seurat_obj)[gene_filter]
  nbin = nbin
  mean_order = order(mvp.mean)
  feature_genes = list()
  gene.bin = floor(gene.n / nbin)
  # gene.bin = seq.int(from = 0, to = 1, length.out = nbin)[2] * gene.n
  for(i in 1:nbin) {
    if(i != nbin) {
      range_in = (gene.bin * (i - 1) + 1) : (gene.bin * i)
    } else {
      range_in = (gene.bin * (i - 1) + 1) : gene.n
    }
    range_in = mean_order[range_in]
    gene_select = gene.list[range_in]
    test_x = mvp.mean[range_in]
    test_y = mvp.dispersion[range_in]
    z_zcore = (test_y - mean(test_y)) / sd(test_y)
    outlier = z_zcore >= 1
    feature_genes = append(feature_genes, gene_select[outlier])
    data_ = data.frame(x = test_x[outlier], y = test_y[outlier])
    figures = figures + geom_point(data=data_, aes(x=x, y=y), size=0.5, color="red")
  }
  return(list(
    figures = figures,
    feature.genes = as.character(feature_genes)
  ))
}
