Rcpp::sourceCpp(paste0(getwd(), "/UCRSI.cpp"))
UCRSI = function(counts, minGene = 4, n_neighbors = 10, outputDistance = F, 
                 alpha_start = 0, alpha_end = 1, alpha_step = 0.1, nJobs = 31, infValue = 0x7fffffff) {
  return(calcDistance(counts, minGene, n_neighbors, outputDistance, alpha_start, 
                      alpha_end, alpha_step, nJobs, infValue))
}

