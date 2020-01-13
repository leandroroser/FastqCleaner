
#' Filter sequences of a FASTQ file by length 
#' @param input  \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param rm.min Threshold value for the minimun number of bases 
#' @param rm.max Threshold value for the  maximum number of bases
#' @description The program removes from a ShortReadQ object those sequences
#' with a length lower than rm.min or/and higher than rm.max
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} 
#' object
#' 
#' @examples 
#' require('Biostrings')
#' require('ShortRead')
#' 
#' # create  ShortReadQ object width widths between 1 and 100
#'  
#'  set.seed(10)
#' input <- random_length(100, widths = 1:100)
#' 
<<<<<<< HEAD
#' # apply the filter, removing sequences length < 10 or length > 80
=======
#' # apply the filter, removing sequences with  10> length > 80
>>>>>>> upstream/master
#' filtered <- length_filter(input, rm.min = 10, rm.max = 80)
#' 
#' # look at the filtered sequences
#' sread(filtered)
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export
#'
#'


length_filter <- function(input, rm.min = NA, rm.max = NA) {
if (!is.na(rm.min) && !is.na(rm.max)){
    if (rm.min > rm.max) {
        stop("The value of rm.min cannot be smaller than rm.max")
    }
}

if (!is.na(rm.max)) {
    input <- input[width(input) <= rm.max]
}
if (!is.na(rm.min)) {
    input <- input[rm.min <= width(input)]
}

input
}
