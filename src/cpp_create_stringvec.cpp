//' create a StringVector using the size of a given NumericVector
//' @param vec Input NumericVector
//' @param x character used to fill the output vector

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11]]

// [[Rcpp::export]]
StringVector cpp_create_stringvec(NumericVector vec, char x) 
{
  size_t vec_size = vec.size();
  StringVector output(vec_size);
  std::string temporal;
  
  for(size_t i = 0; i < vec_size; ++i)
  {
    temporal = "";
    
    for(size_t j = 0; j < vec[i]; ++j)
      temporal += x;
    
    output[i] = temporal;
  }
  return output;
}
