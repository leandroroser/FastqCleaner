#'
#'  Remove a fixed number of bases of a ShortReadQ object from 3' or 5'
#' @param input  \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param trim3 number of bases to remove from 3'
#' @param trim5 number of bases to remove from 5'
#' @description The program removes a given number of bases from the 3' or 5'
#' extremes of the sequences in a ShortReadQ object
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @examples
#' require('Biostrings')
#' require('ShortRead')
#' 
#' # create 6 sequences of width 20
#'  
#' input <- random_seq(6, 20, seed = 10)
#' 
#' # create qualities of width 20 
#' 
#' input_q <- random_qual(c(30,40), slength = 6, swidth = 20, 
#' seed = 10, encod = 'Sanger')
#' 
#' 
#' # create names
#' input_names <- seq_names(6)
#' 
#' # create ShortReadQ object
#' my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)
#' 
#' # apply the filter 
#' filtered3 <- fixed_filter(my_read, trim5 = 5)
#' 
#' filtered5 <- fixed_filter(my_read, trim3 = 5)
#' 
#' filtered3and5 <- fixed_filter(my_read, trim3 = 10, trim5 = 5)
#' 
#' # watch the trimmed sequences
#' sread(filtered3)
#' sread(filtered5)
#' sread(filtered3and5)
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export
#'
#'


fixed_filter <- function(input, trim3 = NA, trim5 = NA) {

if (!is.na(trim3)) {
    input <- input[!(width(input) <= trim3)]
    input <- narrow(input, start = trim3 + 1)
}

if (!is.na(trim5)) {
    input <- input[!(width(input) <= trim5)]
    input <- narrow(input, end = width(input) - trim5)
}

input
}

