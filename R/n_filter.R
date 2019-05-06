
#' Remove sequences with non-identified bases (Ns) from a ShortReadQ object
#' @param input \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param rm.N Threshold value of N's to remove a sequence from the output 
#' (sequences with number of Ns > threshold are removed)
#' For example, if rm.N is 3, all the sequences 
#' with a number of Ns > 3 (Ns >= 4) will be removed
#' @description This program is a wrapper to
#' \code{\link[ShortRead:srFilter]{nFilter}}. 
#' It removes the sequences with a number of N's above 
#' a threshold value 'rm.N'. 
#' All the sequences with a number of N > rm.N (N >= rm.N) will be removed
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} 
#' object
#' @examples 
#' 
#' require('Biostrings')
#' require('ShortRead')
#' 
#' # create 6 sequences of width 20
#' set.seed(10)
#' input <- random_seq(50, 20)
#' 
#' # inject N's
#' set.seed(10)
#' input <- inject_letter_random(input, how_many_seqs = 1:30, 
#' how_many = 1:10)
#'
#' input <- DNAStringSet(input)
#' 

#'
#' # watch the N's frequency
#' hist(letterFrequency(input, 'N'), breaks = 0:10, 
#' main  = 'Ns Frequency', xlab = '# Ns')
#'
#' # create qualities of width 20
#' set.seed(10) 
#' input_q <- random_qual(50, 20)
#'
#' # create names
#' input_names <- seq_names(50)
#'
#' # create ShortReadQ object
#' my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)
#'
#' # apply the filter 
#' filtered <- n_filter(my_read, rm.N = 3)
#'
#' # watch the filtered sequences
#' sread(filtered)
#' 
#' # watch the N's frequency
#' hist(letterFrequency(sread(filtered), 'N'), 
#' main = 'Ns distribution', xlab = '')
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export


n_filter <- function(input, rm.N) {
    ThisFilter <- nFilter(threshold = rm.N)
    input[ThisFilter(sread(input))]
}
