#define ARMA_64BIT_WORD 1
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
#include <RcppThread.h>
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm> 
  
using namespace arma;
using namespace Rcpp;


std::vector<float> sort_indexes(const std::vector<float> &v) {
  std::vector<float> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&v](float i1, float i2) {return v[i1] < v[i2];});
  return idx;
}

int calcIntersectionSum(int** countsBool, int p, int q, int ngene) {
  int sum = 0;
  for(int i = 0; i < ngene; i++) {
    if(countsBool[p][i] > 0 && countsBool[q][i] > 0) {
      sum += 1;
    }
  }
  return sum;
}
float cosine(float* p, float* q, int size) {
  float dot = 0;
  float length_p = 0;
  float length_q = 0;
  for(int i = 0; i < size; i++) {
    dot += p[i] * q[i];
    length_p += p[i] * p[i];
    length_q += q[i] * q[i];
  }
  return dot / (std::sqrt(length_p) * std::sqrt(length_q));
}
float* calcDiff(float* p, float* q, int ngene, int D_pq) {
  float* result = new float[D_pq];
  int index = 0;
  for(int i = 0; i < ngene; i++) {
    if(p[i] > 0 && q[i] > 0) {
      result[index] = p[i] - q[i];
      index ++;
    }
  }
  if(D_pq != index) {
    RcppThread::Rcout << "calcDiff" << std::endl;
  }
  return result;
}
int L0norm(float* diff, int size) {
  int result = 0;
  for(int i = 0; i < size; i++) {
    if(diff[i] != 0) {
      result += 1;
    }
  }
  return result;
}
float L1norm(float* diff, int size) {
  float result = 0;
  for(int i = 0; i < size; i++) {
    result += std::abs(diff[i]);
  }
  return result;
}
float L2norm(float* diff, int size) {
  float result = 0;
  for(int i = 0; i < size; i++) {
    result += diff[i] * diff[i];
  }
  return std::sqrt(result);
}
float arrayMin(float* array, int size) {
  float result = array[0];
  for(int i = 1; i < size; i++) {
    result = std::min(array[i], result);
  }
  return result;
}
float arrayMax(float* array, int size) {
  float result = 0;
  for(int i = 0; i < size; i++) {
    result = std::max(array[i], result);
  }
  return result;
}
float arrayMean(float* array, int size) {
  float result = 0;
  for(int i = 0; i < size; i++) {
    result += array[i];
  }
  return result / size;
}
void maxMinNormalization(float** matrix, const int ncell, int infValue) {
  float distance_min = infValue, distance_max = 0;
  for(int i = 0; i < ncell - 1; i++) {
    for(int i2 = i + 1; i2 < ncell; i2++) {
      distance_min = std::min(distance_min, matrix[i][i2 - i - 1]);
      if(matrix[i][i2 - i - 1] < infValue) {
        distance_max = std::max(distance_max, matrix[i][i2 - i - 1]);
      }
      
    }
  }
  float diff = distance_max - distance_min;
  for(int i = 0; i < (ncell - 1); i++) {
    for(int i2 = i + 1; i2 < ncell; i2++) {
      if(matrix[i][i2 - i - 1] < infValue) {
        matrix[i][i2 - i - 1] = (matrix[i][i2 - i - 1] - distance_min) / diff;
      }
    }
  }
  
}

float mean(float *x, int len) {
  float sum = 0;
  for (int i = 0; i < len; i++)
    sum += x[i];
    return sum / len;
}

float variance(float *x, int len) {
  float meanOfSample = mean(x, len);
  float sum = 0;
  for (int i = 0; i < len; i++) 
    sum += std::pow(x[i] - meanOfSample, 2);
    return sum / (len - 1);
}
float sd(float *x, int len) {
  return std::sqrt(variance(x, len));
}
float varCoff(float *x, int len) {
  return sd(x, len) / mean(x, len);
}

void adjacencyMerge(float** distance, mat* adjacency_final, int size, int n_neighbors = 10, int nJobs = 31) {
  const float minValue = 10;
  RcppThread::ProgressBar bar(size, 1);
  RcppThread::parallelFor(0, size, [adjacency_final, distance, size, n_neighbors, minValue, &bar] (int i) {
    std::vector<float> copy(size);
    for(int i2 = 0; i2 < size; i2++) {
      if(i == i2) {
        copy[i2] = 0;
      } else if(i > i2) {
        copy[i2] = distance[i2][i - i2 - 1];
      } else {
        copy[i2] = distance[i][i2 - i - 1];
      }
    }
    std::vector<float> copy_sort_index = sort_indexes(copy);
    int tIndex = -1;
    for(int i2 = 1; i2 <= n_neighbors; i2++) {
      tIndex = copy_sort_index[i2];
      if(i > tIndex && distance[tIndex][i - tIndex - 1] < minValue) {
        (*adjacency_final)(i, tIndex) += 1 - distance[tIndex][i - tIndex - 1];
        (*adjacency_final)(tIndex, i) = (*adjacency_final)(i, tIndex);
      } else if(i < tIndex && distance[i][tIndex - i - 1] < minValue) {
        (*adjacency_final)(i, tIndex) += 1 - distance[i][tIndex - i - 1];
        (*adjacency_final)(tIndex, i) = (*adjacency_final)(i, tIndex);
      }
    }
    bar++;
  }, nJobs);
}

// [[Rcpp::export]]
NumericMatrix distanceListToMatrix(List distance) {
  const int ncells = distance.size();
  NumericMatrix re(ncells + 1, ncells + 1);
  for(int i = 0; i < ncells; i++) {
    NumericVector v = distance[i];
    for(int i2 = (i + 1); i2 < (ncells + 1); i2++) {
      re(i, i2) = v[i2 - i - 1];
      re(i2, i) = v[i2 - i - 1];
    }
  }
  return re;
}

// [[Rcpp::export]]
List calcDistance(NumericMatrix counts, int minGene = 4, int n_neighbors = 10, 
                  bool outputDistance = false, float alpha_start = 0, float alpha_end = 1, 
                  float alpha_step = 0.1, int nJobs = 31, int infValue = 0x7fffffff) {
  int ncell = counts.ncol();
  int ngene = counts.nrow();
  
  float** countsMatrix = new float* [ncell];
  for(int i = 0; i < ncell; i++) {
    countsMatrix[i] = new float[ngene];
    for(int i2 = 0; i2 < ngene; i2++) {
      countsMatrix[i][i2] = counts(i2, i);
    }
  }
  
  RcppThread::Rcout << "Generating Gene Boolean expression matrix." << std::endl;
  int** countsBool = new int* [ncell];
  for(int i = 0; i < ncell; i++) {
    countsBool[i] = new int[ngene];
    for(int i2 = 0; i2 < ngene; i2++) {
      countsBool[i][i2] = (counts(i2, i) > 0) ? 1 : 0;
    }
  }
  
  RcppThread::Rcout << "Summary non-zero expression." << std::endl;
  int* nonzeroSum = new int[ncell];
  for(int i = 0; i < ncell; i++) {
    int sum = 0;
    for(int i2 = 0; i2 < ngene; i2++) {
      if(countsBool[i][i2] > 0) {
        sum += 1;
      }
    }
    nonzeroSum[i] = sum;
  }
  
  RcppThread::Rcout << "Calculating distance." << std::endl;
  float** distance1 = new float* [ncell - 1];
  float** distance2 = new float* [ncell - 1];
  float** distanceFinal = new float* [ncell - 1];
  for(int i = 0; i < (ncell - 1); i++) {
    distance1[i] = new float[ncell - i - 1];
    distance2[i] = new float[ncell - i - 1];
    distanceFinal[i] = new float[ncell - i - 1];
  }
  
  RcppThread::ProgressBar bar(ncell - 1, 1);
  RcppThread::parallelFor(0, ncell - 1, [ncell, ngene, infValue, minGene, nonzeroSum, 
                          countsBool, distance1, distance2, countsMatrix, &bar] (int i) {
    for(int i2 = (i + 1); i2 < ncell; i2++) {
     int D_pq = calcIntersectionSum(countsBool, i, i2, ngene);
     
     if(D_pq < minGene) {
       distance1[i][i2 - i - 1] = infValue * 1.0;
       distance2[i][i2 - i - 1] = infValue * 1.0;
       continue;
     }
     float c = nonzeroSum[i] * nonzeroSum[i2] * 1.0 / (D_pq * D_pq);
     distance1[i][i2 - i - 1] = c * (1 - cosine(countsMatrix[i], countsMatrix[i2], ngene));
     float* diff_pq = calcDiff(countsMatrix[i], countsMatrix[i2], ngene, D_pq);
     
     float L2 = L2norm(diff_pq, D_pq);
     if(L2 == 0) {
       distance2[i][i2 - i - 1] = 0;
     } else {
       float L1 = L1norm(diff_pq, D_pq);
       distance2[i][i2 - i - 1] = c * L1 * std::sqrt(D_pq) / L2;
     }
     delete[] diff_pq;
     
    }
    bar++;
    
  }, nJobs);
  
  // Min-max normalization
  RcppThread::Rcout << "Min-max normalization." << std::endl;
  float *** distancePackage = new float**[2];
  distancePackage[0] = distance1;
  distancePackage[1] = distance2;
  RcppThread::parallelFor(0, 2, [distancePackage, ncell, infValue] (int i) {
    maxMinNormalization(distancePackage[i], ncell, infValue);
  }, 2);
  delete[] distancePackage;
  
  List resultsList;
  if(outputDistance) {
    RcppThread::Rcout << "Output distance list." << std::endl;
    List distance1_list = List::create();
    List distance2_list = List::create();
    for(int i = 0; i < (ncell - 1); i++) {
      NumericVector v1(ncell - i - 1);
      NumericVector v2(ncell - i - 1);
      for(int i2 = i + 1; i2 < ncell; i2++) {
        v1[i2 - i - 1] = distance1[i][i2 - i - 1];
        v2[i2 - i - 1] = distance2[i][i2 - i - 1];
      }
      distance1_list.push_back(v1);
      distance2_list.push_back(v2);
    }
    resultsList["distance1"] = distance1_list;
    resultsList["distance2"] = distance2_list;
  }

  RcppThread::Rcout << "Merging adjacency matrix." << std::endl;
  mat adjacency_final = mat(ncell, ncell);
  for(float alpha = alpha_start; alpha < alpha_end + alpha_step; alpha+=alpha_step) {
    for(int i = 0; i < (ncell - 1); i++) {
      for(int i2 = i + 1; i2 < ncell; i2++) {
        if(distance1[i][i2 - i - 1] < infValue) {
          distanceFinal[i][i2 - i - 1] = distance1[i][i2 - i - 1] * (1 - alpha) + distance2[i][i2 - i - 1] * alpha;
        } else {
          distanceFinal[i][i2 - i - 1] = 100;
        }
      }
    }
    adjacencyMerge(distanceFinal, &adjacency_final, ncell, n_neighbors, nJobs);
  }
  resultsList["adjacency"] = sp_mat(adjacency_final);
  
  RcppThread::Rcout << "Free memory." << std::endl;
  for(int i = 0; i < ncell; i++) {
    delete[] countsBool[i];
    delete[] countsMatrix[i];
  }
  for(int i = 0; i < (ncell - 1); i++) {
    delete[] distance1[i];
    delete[] distance2[i];
    delete[] distanceFinal[i];
  }
  delete[] distance1;
  delete[] distance2;
  delete[] distanceFinal;
  delete[] countsBool;
  delete[] countsMatrix;
  delete[] nonzeroSum;
  return resultsList;
}
