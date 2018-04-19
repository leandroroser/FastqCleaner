
#' Filter sequences with low quality in 3' tails
#' @param input \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param rm.3qual quality threshold for 3' tails
#' @param q_format quality format used for the file, 
#' as returned by check_encoding
#' @param check.encod check the encoding of the sequence? This argument 
#' is incompatible with q_format. Default TRUE
#' @param remove_zero remove zero-length sequences?
#' @description The program removes the 3' tails with 
#' nucleotides < a qualiy threshold 
#' in a FASTQ file
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} 
#' object
#' @examples 
#'
#' require('Biostrings')
#' require('ShortRead')
#'
#' # create 6 sequences of width 20
#'  
#' input <- random_seq(6, 20, seed = 10)
#' 
#' # create qualities of width 15 and paste to qualities
#' # of length 5 used for the tails.
#' # for two of the sequences, put low qualities in tails
#' 
#' my_qual <- random_qual(c(30,40), slength = 6, swidth = 15, 
#' seed = 10, encod = 'Sanger')
#' tails <-   random_qual(c(30,40), slength = 6, swidth = 5, 
#' seed = 10, encod = 'Sanger')
#' tails[2:3] <- random_qual(c(3, 20), slength = 2, 
#' swidth = 5, seed = 10, encod = 'Sanger')
#' my_qual <- paste0(my_qual, tails)
#' input_q <- BStringSet(my_qual)

#' # create names
#' input_names <- seq_names(6)
#' 
#' # create ShortReadQ object
#' my_read <- ShortReadQ(sread = input,
#' quality = input_q, id = input_names)
#' 
#' # apply the filter 
#' filtered <- trim3q_filter(my_read, rm.3qual = 28)
#' 
#' # watch the trimmed sequences
#' sread(filtered)

#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export


trim3q_filter <- function(input, rm.3qual, q_format = NULL, check.encod = TRUE, 
                        remove_zero = TRUE) {

if (is.null(q_format) && !check.encod) {
    stop("Additional argument must be supplied: 
        q_format must be not null or check.encod must be TRUE")
}

if (!is.null(q_format) && check.encod) {
    check.encod <- FALSE
}


input <- create_uniform_width(input, type = "fastq")
is_edited <- attributes(input)$editedWidth

seqs <- sread(input)
qual <- quality(quality(input))

myqual_mat <- matrix(utf8ToInt(as.character(unlist(qual))), 
nrow = length(qual), 
byrow = TRUE)

if (is_edited) {
    zero_q <- which(myqual_mat == 125)
}

if (!is.null(q_format)) {
    my_encoding <- q_format
    myqual_mat <- myqual_mat - my_encoding$q
}


if (is_edited) {
    myqual_mat[zero_q] <- NA
}

if (check.encod) {
    my_encoding <- check_encoding(myqual_mat)
    q_value <- my_encoding$q
if (my_encoding$y %in% c(1, 2)) {
    myqual_mat <- myqual_mat - 33
} else {
    myqual_mat <- myqual_mat - 64
}
}


if (my_encoding$y %in% c(1, 2, 3)) {
    myqual_mat[is.na(myqual_mat)] <- 0
} else if (my_encoding$y == 4) {
    myqual_mat[is.na(myqual_mat)] <- 3  
} else if (my_encoding$y == 5) {
    myqual_mat[is.na(myqual_mat)] <- -5 
}

at <- myqual_mat < as.integer(rm.3qual)
letter_subject <- DNAString(paste(rep.int("N", ncol(at)), collapse = ""))


letter <- as(IRanges::Views(letter_subject, start = 1,
    end = rowSums(at)), "DNAStringSet")

injectedseqs <- replaceLetterAt(seqs, at, letter)


adapter <- paste(rep("N", max(width(injectedseqs))), sep = "", collapse = "")
mismatchVector <- c(rep(0, width(adapter)))

trimCoords <- trimLRPatterns(Rpattern = adapter,
    subject = injectedseqs,
    max.Rmismatch = mismatchVector, 
    ranges = TRUE)

out <- narrow(input, start = start(trimCoords), end = end(trimCoords))

if (remove_zero) {
    which.zero <- which(width(out) == 0)
if (length(which.zero) > 0) {
    out <- out[-which.zero]
}
}
out
}

