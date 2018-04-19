
//'Compute the base content in a DNA matrix
//'@param x matrix with nucleotides for each position (column)
//'of the reads (rows)
//'@param relative. Logical. If true (default), the function
//' computes relative frequencies

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11]]

// [[Rcpp::export]]
NumericMatrix cpp_base_content(StringMatrix x, bool relative = true ) {
 
 size_t rows = x.nrow();
 size_t cols = x.ncol();
 
StringMatrix bases(1,5);

bases(0, 0) = "C";
bases(0, 1) = "G";
bases(0, 2) = "A";
bases(0, 3) = "T";
bases(0, 4) = "N";


NumericMatrix out(5, cols);

for(size_t j = 0; j < cols; ++j)
{
  for(size_t i = 0; i < rows; ++i)
  {
    if(x(i,j) == bases(0, 0)) {
      out(0, j) += 1;
    } else if(x(i,j) == bases(0, 1)) {
      out(1, j) += 1;
    } else if(x(i,j) == bases(0, 2)) {
      out(2, j) += 1;
    }  else if(x(i,j) == bases(0, 3)) {
      out(3, j) += 1;
    } else if(x(i,j) == bases(0, 4)) {
      out(4, j) += 1;
    }
  }
  if(relative) {
    double total = 0;
    for(size_t k = 0; k < 5; ++k){
      total += out(k, j);
    }
    for(size_t k = 0; k < 5; ++k){
      out(k, j) = out(k, j) / total;
    }
  }
}

return out;
}
