
#' Remove duplicated sequences in a FASTQ file
#' @param input \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @description This program is a wrapper to  
#' \code{\link[ShortRead:srFilter]{occurrenceFilter}}. 
#' It removes the duplicated sequences of a FASTQ file.
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @examples 
#' 
#' require('Biostrings')
#' require('ShortRead')
#' 
#' set.seed(10)
#' s <- random_seq(10, 10)
#' s <- sample(s, 30, replace = TRUE)
#' q <- random_qual(30, 10)
#' n <- seq_names(30)
#'
#' my_read <- ShortReadQ(sread = s, quality = q, id = n)
#' 
#' # check presence of duplicates
#' isUnique(as.character(sread(my_read)))
#' 
#' # apply the filter
#'filtered <- unique_filter(my_read)
#'
#' isUnique(as.character(sread(filtered)))
#' @export

unique_filter <- function(input) {
    uniquef <- occurrenceFilter(withSread = TRUE)
    input[uniquef(input)] 
}
