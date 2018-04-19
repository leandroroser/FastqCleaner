
#' Is natural number
#' @return Logical
#' @keywords  internal

isNaturalNumber <- function(x) {


if (length(x) != length(unlist(x))) {
    return(FALSE)
}
x <- unlist(x)


if (any(is.na(x))) {
    return(FALSE)
}

if (!is.numeric(x)) {
    return(FALSE)
}

if (!all(x == round(x) & x > 0)) {
    return(FALSE)
}
TRUE
}


#' ASCII to integer
#' @return  Integer
#' @keywords internal

asc2int <- function(x) strtoi(charToRaw(x), 16L)


#' integer to ASCII 
#' @return ASCII character
#' @keywords internal

int2asc <- function(n) rawToChar(as.raw(n))


#' check quality encoding
#' @param x quality values
#' @param custom custom encoding from the following:
#' 
#' 'Sanger' --------> expected range: [0, 40]
#' 
#' 'Illumina1.8' --------> expected range: [0, 41]
#' 
#' 'Illumina1.5' --------> expected range: [0, 40]
#' 
#' 'Illumina1.3' --------> expected range: [3, 40]
#' 
#' 'Solexa' --------> expected range: [-5, 40]
#' 
#' @examples
#' 
#' require(Biostrings)
#' 
#' x <- list(PhredQuality(0:40), SolexaQuality(-5:40), IlluminaQuality(3:40))
#' x <- lapply(x, function(i)utf8ToInt(as.character(i)[1]))
#' lapply(x, check_encoding)
#'  
#' SolexaQuality(0:40)
#' IlluminaQuality(0:40)
#' @return  List with encoding information
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export

check_encoding <- function(x = NULL, custom = NULL) {

if (!is.null(x) && !is.null(custom)) {
    stop("Either only a set of qualities 
    or a custom encoding must be selected")
}


if (!is.null(custom) && !custom %in% seq_len(5)) {
    stop("argument custom must be in the range [1, 5]")
}

outFormats <- list(
    A = list(x = "Sanger", y = 1, q = 33, range = "[33 - 73]"), 
    B = list(x = "Illumina 1.8+", y = 2, q = 33, range = "[33; 74]"), 
    C = list(x = "Illumina 1.5+", y = 3, q = 64, range = "[66; 104]"), 
    D = list(x = "Illumina 1.3+", y = 4, q = 64, range = "[64; 104]"),
    E = list(x = "Solexa", y = 5, q = 64, range = "[59; 104]"))

if (!is.null(x)) {
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)

if (xmin < 33 || xmax > 104) {
    stop(paste0("Quality values corrupt. Found [min = ", 
        xmin, " max=", xmax, "],  
    but a range between [33; 104] was expected\n"))
}

if (xmin >= 33 && xmax <= 73) 
    out <- outFormats$A  
else if (xmin >= 33 && xmax <= 74) 
    out <- outFormats$B  
else if (xmin >= 66 && xmax <= 104) 
    out <- outFormats$C 
else if (xmin >= 64 && xmax <= 104) 
    out <- outFormats$D 
else if (xmin >= 59 && xmax <= 104) 
    out <- outFormats$E 

} else {
if (custom == 1) 
    out <- outFormats$A  
else if (custom == 2) 
    out <- outFormats$B 
else if (custom == 3) 
    out <- outFormats$C  
else if (custom == 4) 
    out <- outFormats$D 
else if (custom == 5) 
    out <- outFormats$E 
}

out

}
