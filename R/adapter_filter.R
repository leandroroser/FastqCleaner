
#' Remove full and partial adapters from a ShortReadQ object
#' @param input  \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @param Lpattern 5' pattern (character or
<<<<<<< HEAD
#'  \code{\link[Biostrings]{DNAString}} object)
#' @param Rpattern 3' pattern (character or 
#' \code{\link[Biostrings:DNAString-class]{DNAString}} object)
=======
#'  \code{\link[Biostrings:DNAString-class]{DNAString}} object)
#' @param Rpattern 3' pattern (character or 
#' \code{\link[Biostrings:DNAString-class]{DNAString}} object)
#' @param method  Method used for trimming. If 'exact' the method is based
#' on the exact match between the possible subsequences of the subject and 
#' adapter(s).
#' If 'er' the metod is based on the error-rate between the subsequences, 
#' allowing mismatches in any place
>>>>>>> upstream/master
#' @param rc.L Reverse complement Lpattern? default FALSE
#' @param rc.R Reverse complement Rpatter? default FALSE
#' @param first trim first right('R') or left ('L') side of sequences
#'  when both Lpattern and Rpattern are passed
<<<<<<< HEAD
#' @param error_rate Error rate (value in the range [0, 1] 
=======
#' @param error_rate Error rate (value in the range [0, 1] used for the 
#' 'er' method.
>>>>>>> upstream/master
#' The error rate is the proportion of mismatches allowed between
#' the adapter and the aligned portion  of the subject.
#' For a given adapter A, the number of allowed mismatches between each 
#' subsequence s of A and the subject is computed as: error_rate * L_s,
#' where L_s is the length of the subsequence s
#' @param with_indels Allow indels? This feature is available only 
#' when the error_rate is not null
#' @param anchored Adapter or partial adapter within sequence
<<<<<<< HEAD
#' (anchored = FALSE, default) or only in 3' and 5' terminals? 
#' (anchored = TRUE)
#' @param fixed Parameter passed to 
#' \code{\link[Biostrings]{trimLRPatterns}} 
#' Default 'subject',  ambiguities in the pattern only are interpreted 
#' as wildcard. See the argument fixed in 
#' \code{\link[Biostrings]{trimLRPatterns}} 
#' @param remove_zero Remove zero-length sequences? Default TRUE
#' @param checks Perform checks? Default TRUE
#' @param min_match_flank Do not trim in flanks of the subject,
#'if a match has min_match_flank of less length. Default 1L 
#'(only trim with >=2 coincidences in a flank match)
#' @param ... additional parameters passed to
#' \code{\link[Biostrings]{trimLRpatterns}}
#' @return  Edited \code{\link[Biostrings:DNAString-class]{DNAString}} or 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' @description This program can remove adapters and partial 
#' adapters from 3' and 5', using the functions
#' \code{\link[Biostrings]{trimLRPatterns}} 
#'  The program extends the methodology of
#' the  \code{\link[Biostrings]{trimLRPatterns}} function of \pkg{Biostrings},
#' being also capable of removing adapters present within reads and with other
#' additional otpions 
#' (e.g., threshold of minimum number of bases for trimming). 
#' For a given position in the read, the two Biostrings functions return TRUE 
#' when a match is present between a substring of the read and the adapter.
#' As \code{\link[Biostrings]{trimLRPatterns}} , adapter_filter also selects
#  the longest subsequence that includes the aligned
#' region and goes up to the end of the sequence in the corresponding flank 
#' as the best match. The default error rate is 0.2.
#' If several valid matches are found, the function removes the
#' largest subsequence. Adapters can be anchored or not.  
#' When indels are allowed, the second method uses the 'edit distance' between 
#' the subsequences and the adapter
=======
#' (anchored = FALSE, default) or only in 3' and 5' terminals? (anchored = TRUE)
#' @param fixed Parameter passed to 
#' \code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} 
#' or \code{\link[Biostrings:lowlevel-matching]{isMatchingEndingAt}}. 
#' Default 'subject', where only ambiguities in the pattern are
#'  interpreted as wildcard
#' @param remove_zero Remove zero-length sequences? Default TRUE
#' @param checks Perform checks? Default TRUE
#' @param ... additional parameters passed to  
#' \code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} 
#' or \code{\link[Biostrings:lowlevel-matching]{isMatchingEndingAt}}.
#' @param min_match_flank When a match is found between a subsequence
#' of the subject and the adapter in the corresponding flank,  
#' which would be the minimum length of the overlapping region 
#' (threshold) used for trimming?  Default is 1L (trim when >= 2 
#' base(s) match).
#' @description This program can remove adapters and partial 
#' adapters from 3' and 5', using the functions
#'  \code{\link[Biostrings:lowlevel-matching]{isMatchingEndingAt}} and 
#'  \code{\link[Biostrings:lowlevel-matching]{isMatchingStartingAt}} 
#' of \pkg{Biostrings}.  The program extends the methodology of
#' the \code{\link[Biostrings]{trimLRPatterns}} function of \pkg{Biostrings},
#' being also capable of removing adapters present within reads. 
#' For a given position in the read, the two Biostrings functions return TRUE 
#' when a match is present between a substring of the read and the adapter.
#' As \code{\link[Biostrings]{trimLRPatterns}}, adapter_filter also selects
#  the longest subsequence that includes the aligned
#' region and goes up to the end of the sequence in the corresponding flank 
#' as the best match.
#' If several valid matches are found, the function removes the
#' largest subsequence. Adapters can be anchored or not.  
#' Two methods are available:  one based on the exact matching of the adapter
#' and the reads, and other in an error rate. When indels are allowed, 
#' the second method uses the 'edit distance' between the subsequences 
#' and the adapter
>>>>>>> upstream/master
#' 
#' @examples
#' 
#' require('Biostrings')
#' require('ShortRead')
#' 
#' # create 6 sequences of width 43
#' set.seed(10)
#' input <- random_seq(6, 43)
#' 
<<<<<<< HEAD
#' # add adapter in 3' 
#' adapter <- "ATCGACT"
#' 
#' input <- paste0(input, as.character(DNAString(adapter)))
=======
#' # add adapter in 3' reverse complemented. In read 1, 
#' # it will appear the 5' adapter of read 2 reverse complemented.
#' adapter <- 'ATCGACT'
#' 
#' input <- paste0(input, as.character(reverseComplement(DNAString(adapter))))
>>>>>>> upstream/master
#' input <- DNAStringSet(input)
#' 
#' # create qualities of width 50
#' set.seed(10)
#' input_q <- random_qual(c(30,40), slength = 6, swidth = 50, 
#' encod = 'Sanger')
#'
#' # create names
#' input_names <- seq_names(length(input))
#'
#' # create ShortReadQ object
#' my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)
#'
#' # trim adapter
<<<<<<< HEAD
#' filtered <- adapter_filter(my_read, Rpattern = adapter)
=======
#' filtered <- adapter_filter(my_read, Rpattern = adapter, rc.R = TRUE)
>>>>>>> upstream/master
#'
#' # look at the filtered sequences
#' sread(filtered)
#' 
<<<<<<< HEAD
=======
#' # adapter in the second strand of paired-end reads is reverse-complemented,
#' # with adapter in the end of sequence
#' adapterR <- as.character(reverseComplement(DNAString('ATCGACT')))
#' adapterR <- DNAString(adapterR)
#' inputR <- reverseComplement(input)
#' 
#' 
#' # create qualities of width 50
#' set.seed(10)
#' inputqR <- random_qual(c(30,40), slength = 6, swidth = 50,
#' encod = 'Sanger')
#' 
#' my_readR <- ShortReadQ(sread = inputR, quality = inputqR, id = input_names)
#'
#' # trim adapter
#' filteredR <- adapter_filter(my_readR, Rpattern = adapterR)
#'
#' # look at the filtered sequences
#' sread(filteredR)
>>>>>>> upstream/master
#' 
#' @return  Filtered \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} 
#' object
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export
#'

adapter_filter <- function(input, Lpattern = "", 
<<<<<<< HEAD
    Rpattern = "", rc.L = FALSE, rc.R = FALSE, 
    first = c("R", "L"), with_indels = FALSE, 
    error_rate = 0.2, anchored = TRUE, 
    fixed = "subject", remove_zero = TRUE, 
    checks = TRUE, min_match_flank = 3L, ...) {
    
=======
    Rpattern = "", method = c("exact", "er"), 
    rc.L = FALSE, rc.R = FALSE, first = c("R", "L"), with_indels = FALSE, 
    error_rate = 0.2, anchored = TRUE, 
    fixed = "subject", remove_zero = TRUE, checks = TRUE, 
    min_match_flank = 1L, ...) {
    
    method <- match.arg(method)
>>>>>>> upstream/master
    first <- match.arg(first)
    
    if (is.character(Lpattern)) {
        Lpattern <- sub(" ", "", Lpattern)
        Lpattern <- DNAString(Lpattern)
    }
    
    if (is.character(Rpattern)) {
        Rpattern <- sub(" ", "", Rpattern)
        Rpattern <- DNAString(Rpattern)
    }
    
    if (rc.L) {
        Lpattern <- reverseComplement(Lpattern)
    }
    
    if (rc.R) {
        Rpattern <- reverseComplement(Rpattern)
    }
    

    
    if (length(Lpattern) > 0 && length(Rpattern) > 0) {
        
        if (first == "L") {
            
            out <- cutLseq(subject = sread(input),
                Lpattern = Lpattern, fixed = fixed, 
<<<<<<< HEAD
                error_rate = error_rate, anchored = anchored,
                checks = checks,  min_match_flank = min_match_flank, 
                ranges = TRUE, ...)
            
            out <- cutRseq(subject = out, Rpattern = Rpattern, 
                fixed = fixed,  
=======
                method = method, error_rate = error_rate,
                anchored = anchored, checks = checks, 
                min_match_flank = min_match_flank, ranges = TRUE, ...)
            
            out <- cutRseq(subject = out, Rpattern = Rpattern, 
                fixed = fixed, method = method, 
>>>>>>> upstream/master
                error_rate = error_rate, anchored = anchored, 
                checks = checks, min_match_flank = min_match_flank, 
                ranges = TRUE, ...)
        } else {
            
            out <- cutRseq(subject = sread(input), 
                Rpattern = Rpattern, fixed = fixed, 
<<<<<<< HEAD
                error_rate = error_rate, anchored = anchored, 
                checks = checks, min_match_flank = min_match_flank,
                ranges = TRUE, ...)
            
            out <- cutLseq(subject = out, Lpattern = Lpattern,
                fixed = fixed,  error_rate = error_rate, anchored = anchored, 
=======
                method = method, error_rate = error_rate, 
                anchored = anchored, checks = checks, 
                min_match_flank = min_match_flank, ranges = TRUE, ...)
            
            out <- cutLseq(subject = out, Lpattern = Lpattern,
                fixed = fixed, method = method, 
                error_rate = error_rate, anchored = anchored, 
>>>>>>> upstream/master
                checks = checks, min_match_flank = min_match_flank, 
                ranges = TRUE, ...)
        }
        
        
    } else if (length(Lpattern) == 0 && length(Rpattern) > 0) {

        
        out <- cutRseq(subject = sread(input), 
            Rpattern = Rpattern, fixed = fixed, 
<<<<<<< HEAD
            error_rate = error_rate, anchored = anchored, 
            checks = checks, min_match_flank = min_match_flank,
            ranges = TRUE, ...)
=======
            method = method, error_rate = error_rate, 
            anchored = anchored, checks = checks, 
            min_match_flank = min_match_flank, ranges = TRUE, ...)
>>>>>>> upstream/master
        
    } else if (length(Lpattern) > 0 && length(Rpattern) == 0) {

        out <- cutLseq(subject = sread(input), 
            Lpattern = Lpattern, fixed = fixed, 
<<<<<<< HEAD
            error_rate = error_rate, anchored = anchored, 
            checks = checks, min_match_flank = min_match_flank, 
            ranges = TRUE, ...)
=======
            method = method, error_rate = error_rate, 
            anchored = anchored, checks = checks, 
            min_match_flank = min_match_flank, ranges = TRUE, ...)
>>>>>>> upstream/master
        

    } else if (length(Lpattern) == 0 && length(Rpattern) == 0) 
        {
            
            return(input)
        }  
    
    
    
    out <- narrow(input, start = start(out), end = end(out))
    
    
    if (remove_zero) {
        which.zero <- which(width(out) == 0)
        if (length(which.zero) > 0) {
            out <- out[-which.zero]
        }
    }
    
    out
<<<<<<< HEAD
=======
    
>>>>>>> upstream/master
}
