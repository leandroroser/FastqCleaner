# Auxiliar commands for documentation of the package with roxygen2 in Rstudio

<<<<<<< HEAD
#' @useDynLib FastqCleaner
#' @import  Rcpp
# @importFrom shiny renderDataTable
=======
<<<<<<< HEAD
#' @useDynLib FastqCleaner
#' @import  Rcpp
# @importFrom shiny renderDataTable
=======
#'@useDynLib FastqCleaner
#' @import  Rcpp
#' @import shiny except renderDataTable
>>>>>>> upstream/master
>>>>>>> upstream/master
#' @importFrom graphics hist
#' @importFrom methods as is
#' @importFrom S4Vectors isConstant
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom stats quantile
#' @importFrom shinyBS bsPopover bsModal bsAlert
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> upstream/master
#' @importFrom Biostrings DNAString reverseComplement dinucleotideFrequency PhredQuality quality trimLRPatterns oligonucleotideFrequency DNAStringSet BStringSet replaceLetterAt
#' @importFrom ShortRead sread id ShortReadQ FastqSampler yield nFilter
#' @importFrom IRanges narrow start end width
#' @import htmltools
<<<<<<< HEAD


utils::globalVariables(c("DNAString", "reverseComplement", "sread", "narrow", 
                          "start", "end", "width", "dinucleotideFrequency", "PhredQuality", 
                          "quality", "id", "ShortReadQ", " width", "trimLRPatterns", "FastqSampler", 
                          "yield", "nFilter", "oligonucleotideFrequency", "DNAStringSet", 
                          "BStringSet", "replaceLetterAt", "occurrenceFilter", "method"))

=======


utils::globalVariables(c("DNAString", "reverseComplement", "sread", "narrow", 
                          "start", "end", "width", "dinucleotideFrequency", "PhredQuality", 
                          "quality", "id", "ShortReadQ", " width", "trimLRPatterns", "FastqSampler", 
                          "yield", "nFilter", "oligonucleotideFrequency", "DNAStringSet", 
                          "BStringSet", "replaceLetterAt", "occurrenceFilter", "method"))

=======


>>>>>>> upstream/master
>>>>>>> upstream/master
NULL