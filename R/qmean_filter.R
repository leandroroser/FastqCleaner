
#' Filter sequences by their average quality
#' @param input  \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param minq Quality threshold
#' @param q_format Quality format used for the file, as returned 
#' by check.encoding
#' @param check.encod Check the encoding of the sequence? This argument 
#' is incompatible with q_format
#' @description The program removes the sequences with a quality 
#' lower the 'minq' threshold 
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} 
#' object
#' @examples 
#' 
#' require(ShortRead)
#' 
#' set.seed(10)
#' # create 30 sequences of width 20
#' input <- random_seq(30, 20)
#' 
#' # create qualities of width 20 
#' ## high quality (15 sequences)
#' set.seed(10)
#' my_qual <- random_qual(c(30,40), slength = 15, swidth = 20,
#'                        encod = 'Sanger')
#' ## low quality (15 sequences)
#' set.seed(10)
#' my_qual_2 <-   random_qual(c(5,30), slength = 20,  
#'                            encod = 'Sanger')
#' 
#' # concatenate vectors
#' input_q<- c(my_qual, my_qual_2)
#' 
#' # create names
#' input_names <- seq_names(30)
#' 
#' # create ShortReadQ object
#' my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)
#' 
#' # watch the average qualities
#' alphabetScore(my_read) / width(my_read)
#'
#' # apply the filter
#' filtered <- qmean_filter(my_read, minq = 30)
#'
#' # watch the average qualities
#' alphabetScore(my_read) / width(my_read)
#'
#'
#' # watch the filtered sequences
#' sread(filtered)
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export


qmean_filter <- function(input, minq, q_format = NULL, check.encod = TRUE) {

if (is.null(q_format) && !check.encod) {
    stop("Additional argument must be supplied:
    q_format must be not null or check.encod must be TRUE")
}

if (!is.null(q_format) && check.encod) {
    check.encod <- FALSE
}

input <- create_uniform_width(input, type = "fastq")

is_edited <- attributes(input)$editedWidth

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


input[rowMeans(myqual_mat, na.rm = TRUE) >= minq]
}
