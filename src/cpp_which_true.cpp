
//' Select the correct match position for matching operations with
//' multiple strings 
//' @description The program takes as input a logical matrix that 
//' indicates with TRUE or FALSE if a match is present in a given position
//'  or not, respectively. Each column is a sequence, each row a cycle
//'  for NGS experiments
//' @param m LogicalMatrix with TRUE and FALSE values
//' @param width length of each sequence. This allow to discard fake matches
//' @param origin If origin = "start", select the first match. 
//' If origin = "end", selects the last match.
//' If origin  = "end_inverted" selects the first match starting 
//' from n_i positions, where n is the  length of the sequence represented 
//' in the column i.  This case works when the last match is 
//' seeked in reversed sequences.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11]]


// [[Rcpp::export]]
IntegerVector cpp_which_true(LogicalMatrix m, IntegerVector width, String origin)
{
  size_t mcol = m.ncol();
  IntegerVector out(mcol);
  for(size_t i = 0; i < mcol; ++i) {
    out[i] = -1;
  }
  
  if(origin == "start") {
    for(size_t i = 0; i < mcol; ++i)
    {
      for(int j = 0; j < width[i]; ++j)
      {
        if(m(j, i) == true)
        {
          out[i] = j + 1;
          break;
        } 
      }
    }
  }
  
  if(origin == "end"  || origin == "end_inverted") {
  
  for(size_t i = 0; i < mcol; ++i)
  {
    for(int j = 0; j < width[i]; ++j)
    {
      if(m(j, i) == true)
      {
        out[i]  = j + 1;
      }
   }
  }
  }
  
  if(origin == "end_inverted") {
     for(size_t i = 0; i < mcol; ++i) {
       if(out[i] != -1) out[i] = width[i] - out[i] + 1;
     }
  }
  
  return(out);
}
