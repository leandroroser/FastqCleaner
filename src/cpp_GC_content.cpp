//' Compute CG content
//' @param x StringMatrix with bases, of reads x cycles


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11]]

// [[Rcpp::export]]
NumericMatrix cpp_GC_content(StringMatrix x) {
 
 size_t rows = x.nrow();
 size_t cols = x.ncol();
 
StringMatrix bases(1,4);

bases(0, 0) = "C";
bases(0, 1) = "G";
bases(0, 2) = "A";
bases(0, 3) = "T";

NumericMatrix out(rows, 4);

for(size_t i = 0; i < rows; ++i)
{
  for(size_t j = 0; j < cols; ++j)
  {
    if(x(i,j) == bases(0, 0)) {
      out(i, 0) += 1;
    } else if(x(i,j) == bases(0, 1)) {
      out(i, 1) += 1;
    } else if(x(i,j) == bases(0, 2)) {
      out(i, 2) += 1;
    }  else if(x(i,j) == bases(0, 3)) {
      out(i, 3) += 1;
    }
  }
}

return out;
}
