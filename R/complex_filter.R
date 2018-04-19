#' Remove sequences with low complexity
#' @param input  \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param threshold computed H/ reference H treshold
#' @param referenceEntropy reference entropy. By default, 
#' the program uses a value of 3.908, that corresponds 
#' to the entropy of the human genome in bits
#' @description The program removes low complexity sequences, computing the 
#' entropy with the frequency of dinucleotides.
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @examples 
#' 
#' require('Biostrings')
#' require('ShortRead')
#' 
#' # create  sequences of different width
#'  
#'input <- lapply(c(0, 6, 10, 16, 20, 26, 30, 36, 40), 
#'                function(x) random_seq(1, x, seed = 10))
#'
#'
#' # create repetitive 'CG' sequences with length adequante 
#' # for a total length:
#' # input +  CG = 40
#' 
#'CG <- lapply(c(20, 17, 15, 12, 10, 7, 5, 2, 0), 
#'             function(x) paste(rep('CG', x), collapse = ''))
#'
#'
#' # concatenate input and CG
#'input  <- mapply('paste', input, CG, sep = '')
#'input <- DNAStringSet(input)
#'
#' # plot relative entropy (E, Shannon 1948)
#' 
#' freq <- dinucleotideFrequency(input)
#' freq  <- freq /rowSums(freq)
#' H <- -rowSums(freq  * log2(freq), na.rm = TRUE)
#' H_max <- 3.908135  # max entropy
#' plot(H/H_max, type='b', xlab = 'Sequence', ylab= 'E')
#' 
#' 
#' # create qualities of widths 40
#' 
#' input_q <- random_qual(c(30,40), slength = 9, swidth = 40, 
#'                        seed = 10, encod = 'Sanger')
#' 
#' # create names
#' input_names <- seq_names(9)

#' 
#' # create ShortReadQ object
#' my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)
#' 
#' # apply the filter, 
#' filtered <- complex_filter(my_read)
#' 
#' # watch the filtered sequences
#' sread(filtered)
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export
#'

complex_filter <- function(input, threshold = 0.5, 
                                referenceEntropy = 3.908135) {
diNucFreq <- dinucleotideFrequency(sread(input))
if (is.null(dim(diNucFreq))) {
    diNucFreq <- diNucFreq/sum(diNucFreq)
    H <- -sum(diNucFreq * log2(diNucFreq), na.rm = TRUE)
} else {
    diNucFreq <- diNucFreq/rowSums(diNucFreq)
    H <- -rowSums(diNucFreq * log2(diNucFreq), na.rm = TRUE)
}
input[H/referenceEntropy >= threshold]
}
