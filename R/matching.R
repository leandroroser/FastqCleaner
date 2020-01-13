
#' Remove left and right full and partial patterns
#' @param subject \code{\link[Biostrings:DNAString-class]{DNAString}} or 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' @param Lpattern 5' pattern,
#'  \code{\link[Biostrings:DNAString-class]{DNAString}} object
#' @param Rpattern 3' pattern,
#'  \code{\link[Biostrings:DNAString-class]{DNAString}} object
#' @param error_rate Error rate (value in [0, 1]). 
#' The error rate is the proportion of mismatches allowed between
#' the adapter and the aligned portion of the subject.
#' For a given adapter A, the number of allowed mismatches between each 
#' subsequence s of A and the subject is computed as: error_rate * L_s,
#' where L_s is the length of the subsequence s.
#' @param with.indels Allow indels? 
#' @param anchored Can the adapter or partial adapter be within 
#' the sequence? (anchored = FALSE)
#' or only in the terminal regions of the sequence? (anchored = TRUE).
#' Default TRUE (trim only flanking regions)
#' @param fixed Parameter passed to 
#' \code{\link[Biostrings]{trimLRpatterns}}
#' Default 'subject',  ambiguities in the pattern only are interpreted 
#' as wildcard. See the argument fixed in 
#' \code{\link[Biostrings]{trimLRpatterns}}
#' @param ranges Return ranges? Default FALSE
#' @param checks Perform internal checks? Default TRUE
#' @param min_match_flank Do not trim in flanks of the subject,
#'  if a match has min_match_flank of less length. Default 1L 
#'  (only trim with >=2 coincidences in a flank match)
#' @param ... additional parameters passed to
#'  \code{\link[Biostrings]{trimLRpatterns]}}
#' @return  Edited \code{\link[Biostrings:DNAString-class]{DNAString}} or 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' @description This set of programs are internal, 
#' and the function adapter_filter is recommended for trimming. 
#' The programs can remove adapters and partial 
#' adapters from 3' and 5'. The adapters can be anchored or not. 
#' When indels are allowed, the error rate consists in the edit distance. 
#' IUPAC simbols are allowed. The methods use the 
#' \code{\link[Biostrings]{trimLRpatterns}} function
#'  of the \pkg{Biostrings} package, with some additions 
#'  to take into account e.g., partial adaptors.
#'  IUPAC symbols are allowed in all the cases. The present function 
#' also removes partial adapters, without the need of additional steps
#'  (for example, creating a padded adapter with 'Ns', etc). 
#' A similar result to the output of \code{\link[Biostrings]{trimLRPatterns}}
#' can be obtained with the option anchored = TRUE.
#' When several matches are found, the function removes the subsequence 
#' that starts in the first match when cutRseq is used, or ends 
#' in the last match when cutLseq is used.
#' 
#' @examples
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
#' cutLseq(subject, Lpattern,  error_rate = 0.2)
#' cutLseq(subject, Lpattern,  error_rate = 0.2, 
#' with.indels = TRUE)
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @rdname matching
#' @keywords internal


cutRseq <- function(subject, Rpattern, 
                    with.indels = FALSE,
                    fixed = "subject", 
                    error_rate = 0.2, 
                    anchored = TRUE, 
                    ranges = FALSE,
                    checks = TRUE, 
                    min_match_flank = 2L, 
                    ...) {
    
   Rpattern <- DNAString(Rpattern)
  
   if (error_rate > 1 || error_rate < 0) {
    stop("error_rate must be a number between 0 and 1")
   }

   if (checks) {
    
       if (!is(Rpattern, "DNAString")) {
           stop("Rpattern must be a character string or a DNAString object")
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
  
   if(error_rate > 0) {
       flank_seq <- as.integer(seq_len(p) * error_rate)
   } else {
       flank_seq <- rep(0, length(seq_len(p)))
   }
  
  
   if (min_match_flank >= 1L) {
       if (p > min_match_flank) {
           flank_seq[seq_len(min_match_flank)] <- -1
       } else {
           return(subject)
       }
   }
 
   if(!anchored) {
       Rpattern <- as.character(Rpattern)
       maxlen <-  max(width(subject)) - nchar(Rpattern)
       if(maxlen > 0) {
           Rpattern <- paste0(Rpattern,  paste(rep("N",maxlen), collapse = ""))
       }
       Rpattern <- DNAString(Rpattern)
       flank_seq <- c(flank_seq, rep(0,maxlen))
   }
  
   out <- trimLRPatterns(Rpattern = Rpattern,
                       subject = subject, 
                       max.Rmismatch = flank_seq, 
                       with.Rindels = with.indels, 
                       Rfixed = fixed, 
                       ...)
  
  
   if (ranges) {
       out <- IRanges::IRanges(start = rep(1, length(out)), end = width(out))
   } 
  
   out
}

#' Remove left and right full and partial patterns
#' @rdname matching

cutLseq <- function(subject, Lpattern,
                    with.indels = FALSE, 
                    fixed = "subject", 
                    error_rate = 0.2, 
                    anchored = TRUE, 
                    ranges = FALSE, 
                    min_match_flank = 3L, 
                    checks = TRUE, ...) {


   Lpattern <- DNAString(Lpattern)

   if (checks) {
       if (!is(Lpattern, "DNAString")) {
           stop("Rpattern must be a character string or a DNAString object")
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


   if(error_rate > 0) {
       flank_seq <- as.integer(seq_len(p) * error_rate)
   } else {
       flank_seq <- rep(0, length(seq_len(p)))
   }


   if (min_match_flank >= 1L) {
       if (p > min_match_flank) {
           flank_seq[seq_len(min_match_flank)] <- -1
       } else {
           return(subject)
   }
   }


   if(!anchored) {
       Lpattern <- as.character(Lpattern)
       maxlen <-  max(width(subject)) - nchar(Lpattern)
       if(maxlen > 0) {
           Lpattern <- paste0(paste(rep("N",maxlen), collapse = ""), Lpattern)
       }
       Lpattern <- DNAString(Lpattern)
       flank_seq <- c(flank_seq, rep(0,maxlen))
  }


   out <- trimLRPatterns(Lpattern = Lpattern,
                       subject = subject, 
                       max.Lmismatch = flank_seq, 
                       with.Lindels = with.indels, 
                       Lfixed = fixed, 
                       ...)


   if (ranges) {
       out <- IRanges::IRanges(start = rep(1, length(out)), end = width(out))
   } 

   out
}
