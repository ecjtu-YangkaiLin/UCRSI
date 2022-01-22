UCRSI
===
A clustering method Unifying cell-type recognition and subtypes identification for tumor heterogeneity analysis
___

# 1. Install R dependencies package
    install.packages("devtools")
    install.packages("Rcpp")
    install.packages("RcppArmadillo")
    install.packages("Seurat")
    install.packages("ggplot2")
    install.packages("Matrix")
    install.packages("reticulate")
    install.packages("umap")
    devtools::install_github("tnagler/RcppThread")

# 2. Python dependency
UCRSI relies on two Python packages, **scikit-learn** and **pyamg**, to complete the clustering process. In R, we use the **reticulate** package to call it. It is recommended to create a separate virtual environment using **virtualenv**. After configuring the environment, you need to change the **python_env** variable in the **main.R** file to the path of the environment.

# Example
    See main.R
