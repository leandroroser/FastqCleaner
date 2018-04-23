#' Remove a set of sequences
#' @description Removes a set of sequences 
#' @param  input \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param rm.seq Ccharacter vector with sequences to remove
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} 
#' object
#' @examples 
#' 
#' require(ShortRead)
#' 
#' input <- random_length(30, 3:7, seed = 10)
#' rm.seq  = c('TGGTC', 'CGGT', 'GTTCT', 'ATA')
#' 
#' # verify that some sequences match
#' match_before <- unlist(lapply(rm.seq,
#'  function(x) grep(x, as.character(sread(input)))))
#' 
#' filtered <- seq_filter(input,rm.seq =  rm.seq)
#' 
#' # verify that matching sequences were removed
#' match_after <- unlist(lapply(rm.seq, 
#' function(x) grep(x, as.character(sread(filtered)))))
#' 
#' @export

seq_filter <- function(input, rm.seq) {
reads <- sread(input)
matchseq <- lapply(rm.seq, function(x) {
    Biostrings::which.isMatchingStartingAt(x, subject = reads)
})
matchseq <- unlist(lapply(matchseq, function(x) which(x == 1)))
if (length(matchseq) > 0) {
    return(input[-matchseq])
}
input
}
