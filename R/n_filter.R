
#' Remove sequences with non identified bases (N) from a ShortReadQ object
#' 
#' @param input \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param rm.N threshold of N's to remove a sequence in the output 
#' (sequences with N > threshold are removed)
#' file. For example, if rm.N is 3, all the sequences 
#' with a number of N > 3 (N >= 4) will be removed
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
#' input <- random_seq(50, 20, seed = 10)
#' 
#' # inject N's
#' 
#' input <- inject_letter_random(input, how_many_seqs = 1:30, 
#' how_many = 1:10, seed = 10)
#'
#' input <- DNAStringSet(input)
#' 

#'
#' # watch the N's frequency
#' hist(letterFrequency(input, 'N'), breaks = 0:10, 
#' main  = 'Ns Frequency', xlab = '# Ns')
#'
#' # create qualities of width 20 
#' input_q <- random_qual(50, 20, seed = 10)
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

