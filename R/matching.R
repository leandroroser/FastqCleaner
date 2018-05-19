
#' Remove left and right full and partial patterns
#' @param subject \code{\link[Biostrings:DNAString-class]{DNAString}} or 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' @param Lpattern 5' pattern,
#'  \code{\link[Biostrings:DNAString-class]{DNAString}} object
#' @param Rpattern 3' pattern,
#'  \code{\link[Biostrings:DNAString-class]{DNAString}} object
#' @param method  Method used for trimming. If 'exact' the metod is based
#' on the exact matching of the posible subsequences of the subject and 
#' the adapters. If 'er' the metod is based on the eror-rate between 
#' the subsequences, allowing mismatches.
#' @param error_rate Error rate (value in [0, 1] used for 'er' method). 
#' The error rate is the proportion of mismatches allowed between
#' the adapter and the aligned portion of the subject.
#' For a given adapter A, the number of allowed mismatches between each 
#' subsequence s of A and the subject is computed as: error_rate * L_s,
#' where L_s is the length of the subsequence s.
#' @param with.indels Allow indels? This feature is only available
#' for er method.
#' @param anchored Can the adapter or partial adapter be within 
#' the sequence? (anchored = FALSE)
#' or only in the terminal regions of the sequence? (anchored = TRUE).
#' Default TRUE (trim only flanking regions)
#' @param fixed Parameter passed to 
#' code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} and
#' \code{\link[Biostrings:lowlevel-matching]{isMatchingEndingAt}}  
#' Default 'subject',  ambiguities in the pattern only are interpreted 
#' as wildcard. See the argument fixed in 
#' code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} and
#' \code{\link[Biostrings:lowlevel-matching]{isMatchingEndingAt}} 
#' @param ranges Return ranges? Default FALSE
#' @param checks Perform internal checks? Default TRUE
#' @param min_match_flank Do not trim in flanks of the subject,
#'  if a match has min_match_flank of less length. Default 1L 
#'  (only trim with >=2 coincidences in a flank match)
#' @param ... additional parameters passed to
#'  \code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} and
#' \code{\link[Biostrings:lowlevel-matching]{isMatchingEndingAt}} 
#' @return  Edited \code{\link[Biostrings:DNAString-class]{DNAString}} or 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' @description This set of programs are internal, 
#' and the function adapter_filter is recommended for trimming. 
#' The programs can remove adapters and partial 
#' adapters from 3' and 5'. The adapters can be anchored or not. 
#' Two methods are available: one based on the exact matching 
#' of the sequence and the adapter, and a second in an error rate. 
#' When indels are allowed, the error rate consists in the edit distance. 
#' IUPAC simbols are allowed. The methods use the 
#' \code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} function
#'  of the \pkg{Biostrings} package to find matches. IUPAC symbols
#' are allowed in all the cases. The present function 
#' also removes partial adapters, without the need of additional steps
#'  (for example, creating a padded adapter with 'Ns', etc). 
#' A similar result to the output of \code{\link[Biostrings]{trimLRPatterns}}
#' can be obtained with the option anchored = TRUE.
#' When several matches are found, the function removes the subsequence 
#' that starts in the first match when cutRseq is used, or ends 
#' in the last match when cutLseq is used.
#' 
#' @examples
#' \dontrun{
#' library(Biostrings)
#' 
#' subject <- DNAStringSet(c('ATCATGCCATCATGAT',
#' 'CATGATATTA', 'TCATG', 'AAAAAA', 'AGGTCATG'))
#' 
#' Lpattern <- Rpattern <- 'TCATG'
#' 
#' cutLseq(subject, Lpattern)
#' cutLseq(subject, Lpattern, ranges = TRUE)
#' cutRseq(subject, Rpattern)
#' 
#' 
#' cutLseq(subject, Lpattern, anchored = FALSE)
#' cutLseq(subject, Lpattern, method = 'er', error_rate = 0.2)
#' cutLseq(subject, Lpattern, method = 'er', error_rate = 0.2, 
#' with.indels = TRUE)
#' }
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @rdname matching
#' @keywords internal


cutRseq <- function(subject, Rpattern, with.indels = FALSE,
                    fixed = "subject", method = c("exact", "er"), 
                    error_rate = 0.2, 
                    anchored = TRUE, ranges = FALSE,
                    checks = TRUE, min_match_flank = 0L, 
                    ...) {

method <- match.arg(method)

Rpattern <- DNAString(Rpattern)

if (error_rate > 1 || error_rate < 0) {
    stop("error_rate must be a number between 0 and 1")
}


if (checks) {

if (!is(Rpattern, "DNAString")) {
        stop("Rpattern must be a character string or a DNAString object")
}

if (method == "exact") {
    if (with.indels) {
        stop("indels allowed only with er method")
    }
}

csub <- class(subject)
if (csub != "DNAStringSet") {
    stop("subject must be a DNAString or DNAStringSet object")
}

if (csub == "DNAString") {
    subject <- as(subject[[1]], "DNAStringSet")
}

}

p <- length(Rpattern)
s_width <- width(subject)
s <- max(width(subject))

if (method == "er") {
    flank_seq <- rev(as.integer(seq_len(p) * error_rate))
} else {
    flank_seq <- seq_len(p) - 1L
}

fill <- ifelse(anchored, -1L, ifelse(method == "exact", 0L, max(flank_seq)))


maxmis <- c(rep(fill, s - p), flank_seq)

if (min_match_flank >= 1L) {
    if (p > min_match_flank) {
        maxmis[(length(maxmis) - min_match_flank + 1):length(maxmis)] <- -1
    } else {
        return(subject)
    }
}


if (!S4Vectors::isConstant(s_width)) {
which.match <- isMatchingEndingAt(pattern = reverse(Rpattern), 
subject = reverse(subject), 
ending.at = seq_len(s), 
max.mismatch = rev(maxmis),
with.indels = with.indels, 
fixed = fixed, ...)


u <- cpp_which_true(which.match, s_width, "end_inverted")

end.pos <- u - 1


end.pos[u == -1] <- s_width[u == -1]

} else {
which.match <- isMatchingStartingAt(pattern = Rpattern,
subject = subject, 
starting.at = seq_len(s),
max.mismatch = maxmis, 
with.indels = with.indels, 
fixed = fixed, ...)



u <- cpp_which_true(which.match, s_width, "start")

end.pos <- u - 1

end.pos[u == -1] <- s_width[u == -1]


end.pos[u == 1] <- 0

}

start.pos <- rep(1L, length(subject))

if (ranges) {
    return(IRanges::IRanges(start = start.pos, end = end.pos))
} else {
    return(subseq(subject, start = start.pos, end = end.pos))
}

}

#' @rdname remove-adapters
#' @export

cutLseq <- function(subject, Lpattern, with.indels = FALSE, 
                    fixed = "subject", method = c("exact", "er"), 
                    error_rate = 0.2, anchored = TRUE, 
                    ranges = FALSE, min_match_flank = 1L, 
                    checks = TRUE, ...) {

method <- match.arg(method)
Lpattern <- DNAString(Lpattern)

if (checks) {
    if (!is(Lpattern, "DNAString")) {
        stop("Rpattern must be a character string or a DNAString object")
    }

if (method == "exact") {
    if (with.indels) {
        stop("indels allowed only with er method")
    }
}

csub <- class(subject)
if (csub != "DNAStringSet") {
    stop("subject must be a DNAStringSet object")
}

if (csub == "DNAString") {
    subject <- as(subject[[1]], "DNAStringSet")
    }
}



p <- length(Lpattern)
s_width <- width(subject)
s <- max(width(subject))

if (method == "er") {
    flank_seq <- as.integer(seq_len(p) * error_rate)
} else {
    flank_seq <- rev(seq_len(p) - 1L)
}

fill <- ifelse(anchored, -1L, ifelse(method == "exact", 0L, max(flank_seq)))


maxmis <- c(flank_seq, rep(fill, s - p))

if (min_match_flank >= 1L) {
    if (p > min_match_flank) {
        maxmis[seq_len(min_match_flank)] <- -1
} else {
    return(subject)
    }
}

which.match <- isMatchingEndingAt(pattern = Lpattern, 
subject = subject, ending.at = seq_len(s), 
max.mismatch = maxmis, with.indels = with.indels, fixed = fixed, ...)



u <- cpp_which_true(which.match, s_width, "end") 

start.pos <- u + 1L
end.pos <- s_width


start.pos[u == -1] <- 1L


start.pos[start.pos > s_width] <- 1L
end.pos[u == s_width] <- 0L

if (ranges) {
    return(IRanges::IRanges(start = start.pos, end = end.pos))
} else {
    return(subseq(subject, start = start.pos, end = end.pos))
    }
}

