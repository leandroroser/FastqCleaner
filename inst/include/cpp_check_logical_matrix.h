
#ifndef CLM_H
#define CLM_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11]]

// check logical condition with element of a matrix and returns a binary object
// NA values are returned as naReturn

// [[Rcpp::export]]
NumericMatrix cpp_check_logical_matrix(NumericMatrix x, size_t rows, size_t cols, int qThreshold, int naReturn) 
  {
  
  NumericMatrix out(rows, cols);
  for(size_t i = 0; i < rows; ++i) 
    {
    for(size_t j = 0; j < cols; ++j)
      {
      if(!(std::isnan(x(i, j)))) 
        {
          if(x(i, j) >= qThreshold) 
             out(i, j) = 1;
           else 
             out(i, j) = 0;
        } else 
       out(i, j) = naReturn;
  // std::nan("");
      }
    }
  
  return out;
  }

#endif
