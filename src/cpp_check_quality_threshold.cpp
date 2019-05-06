

//' check by row if the values of a matrix are above a 
//' threshold of values in a vector returns a matrix of rows by threshold values
//' @param x NumericMatrix to be checked
//' @param thresVector vector with thresholds

#include <Rcpp.h>
using namespace Rcpp;
#include <cpp_check_logical_matrix.h>

// [[Rcpp::plugins(cpp11]]
// [[Rcpp::export]]
NumericMatrix cpp_check_quality_threshold(NumericMatrix x, 
                                            NumericVector thresVector) {

 size_t rows = x.nrow();
 size_t cols = x.ncol();
 size_t thres_length = thresVector.length();
 NumericMatrix q_mat(rows, thres_length);
 NumericMatrix q_temp_mat(rows, cols);
 int stopval = false;
 
 for(size_t k = 0; k < thres_length; ++k) 
 {
   q_temp_mat = cpp_check_logical_matrix(x, rows, cols, thresVector[k], 1);

   for(size_t i = 0;  i < rows; ++ i) 
   {
     for(size_t j = 0; j < cols; ++ j) 
     {
       if(q_temp_mat(i, j) == 0)
       {
         stopval = true;
         break;
       }
     }
     
     if(!stopval)
       q_mat(i, k) = 1;
     
     stopval = false;
     
   }
 }
 
 return q_mat;
}
