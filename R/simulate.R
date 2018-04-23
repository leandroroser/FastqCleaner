
#' Create random sequences
#' @description Create a
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' with random sequences
#' @param slength Number of sequences
#' @param swidth Width of the sequences
#' @param seed Random seed. Default NULL
#' @param nuc Create sequences of DNA (nucleotides = c('A', 'C', 'G', 'T')) or
#'  RNA (nucleotides = c('A, 'C', 'G', 'U'))?. Default: 'DNA'
#' @param prob A vector of four probability values used 
#' to set the frequency of the nucleotides 'A', 'C', 'G', 'T', for DNA,
#' or 'A', 'C', 'G', 'U', for RNA. For example = c(0.25, 0.25, 0.5, 0).
#' Default is = c(0.25, 0.25, 0.25, 0.25) (equiprobability for the 4 bases).
#' If the sum of the probabilities is > 1, the values will be nomalized 
#' to the range [0, 1].
#' @return \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object
#' 
#' @examples 
#' 
#' s1 <- random_seq(slength = 10, swidth = 20, seed = 10)
#' s1
#' 
#' s2 <- random_seq(slength = 10, swidth = 20, seed = 10, 
#' prob = c(0.6, 0.1, 0.3, 0))
#' s2
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export

random_seq <- function(slength, swidth, seed = NULL, nuc = c("DNA", "RNA"), 
                        prob = c(0.25, 
                        0.25, 0.25, 0.25)) {
    
    
    if (!is.null(seed)) 
        set.seed(seed)
    
    nuc <- match.arg(nuc)
    if (length(prob) > 4 || !is.numeric(prob) || any(prob < 0)) {
        stop("prob must be a numeric vector of length 4")
    }
    
    if (swidth < 0 || length(slength) != 1) {
        stop("swidth must be > 1; slength must be of length 1")
    }
    
    if (length(swidth) > 1) {
        if (length(swidth) != slength) 
            stop("length(swidth) and slength must be equal")
        
    } else {
        swidth <- rep(swidth, slength)
    }
    
    if (nuc == "DNA") {
        samplefrom <- c("A", "C", "G", "T")
    } else if (nuc == "RNA") {
        samplefrom <- c("A", "C", "G", "U")
    }
    
    
    prob <- prob/sum(prob)

    
    out <- list()
    for (i in seq_len(slength)) {
        samp <- sample(samplefrom, swidth[i], replace = TRUE, prob = prob)
        out[[i]] <- paste(samp, collapse = "")
    }
    DNAStringSet(do.call("c", out))
}


#' Create random qualities for a given encoding
#' @param qual quality range for the sequences. It must be a range included 
#' in the selected encoding:
#' 
#' 'Sanger' = [0, 40]
#' 
#' 'Illumina1.8' = [0, 41]
#' 
#' 'Illumina1.5' = [0, 40]
#' 
#' 'Illumina1.3' = [3, 40]
#' 
#' 'Solexa' = [-5, 40]
#'
#' example: for a range from 20 to 30 in Sanger encoding,
#'  pass the argument = c(20, 30)
#' @param slength number of sequences
#' @param swidth width of the sequences
#' @param encod sequence encoding
#' @param seed random seed. Default: NULL
#' @param prob a vector of range = range(qual), with probabilities to set 
#' the frequency of each quality value. Default is equiprobability.
#' If the sum of the probabilities is > 1, the values will be nomalized 
#' to the range [0, 1].
#' @description Create a \code{\link[Biostrings:XStringSet-class]{BStringSet}} 
#' object
#' with random qualities
#' @return  \code{\link[Biostrings:XStringSet-class]{BStringSet}} object
#' @examples 
#' 
#' q <- random_qual(30, 20)
#' q
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export

random_qual <- function(slength, swidth, qual = NULL,
    encod = c("Sanger", "Illumina1.8", 
    "Illumina1.5", "Illumina1.3", "Solexa"), seed = NULL, prob = NULL) {
    
    if (!is.null(seed)) 
        set.seed(seed)
    
    
    encod <- match.arg(encod)
    
    
    if (swidth < 0 || length(slength) != 1) {
        stop("swidth must be > 1; slength must be of length 1")
    }
    
    if (length(swidth) > 1) {
        if (length(swidth) != slength) 
            stop("length(swidth) and slength must be equal")
        
    } else {
        swidth <- rep(swidth, slength)
    }
    
    
    
    selection <- switch(encod, Sanger = 1, Illumina1.8 = 2,
        Illumina1.5 = 3, Illumina1.3 = 4, 
        Solexa = 5, stop("Invalid encoding"))
    

    seqrange <- list(0:40, 0:41, 0:40, 3:40, -5:40)
    
    seqrange <- seqrange[[selection]]
    
    if (!is.null(qual)) {
        if (!all(qual %in% seqrange)) {
            stop(paste0("Quality values out of range. 
                The values must be in the range [", 
                min(seqrange), ",", max(seqrange), "]"))
        }
    } else {
        
        qual <- seqrange
    }

    
    my_scores <- seq(min(qual), max(qual)) + c(33, 33, 33, 64, 64)[selection]
    samplefrom <- vapply(my_scores, int2asc, character(1))
    
    lscores <- length(my_scores)
    if (is.null(prob)) {
        prob <- rep(1, lscores)/lscores
    } else {
    if (length(prob) == lscores || !is.numeric(prob) || any(prob < 0)) {
    stop(paste0("prob must be a positive numeric vector of length swidth (", 
        swidth, ")"))
        }
    }

    
    out <- vapply(seq_len(slength), function(x) {
        paste(sample(samplefrom, swidth, replace = TRUE, prob = prob), 
                    collapse = "")
    }, character(1))
    
    out <- list()
    for (i in seq_len(slength)) {
        out[[i]] <- paste(sample(samplefrom, swidth[i], replace = TRUE, 
                                prob = prob), 
            collapse = "")
    }
    BStringSet(do.call("c", out))
}


#' Create fastq/sequences/qualities with uniform width
#' @param input input to edit
#' @param type type of the input: 'fastq' (ShortReadQ),
#'  'sequence' (DNAStringSet), 'quality' (BStringset)
#' @return ShortReadQ object or character vector with sequences 
#' or qualities, with uniform widht (padded with Ns or \})
#' @keywords internal

create_uniform_width <- function(input,
                                type = c("fastq", "sequence", "quality")) {
    
    if (type == "fastq") {
        
        seqs <- sread(input)
        qual <- PhredQuality(quality(quality(input)))
        ids <- id(input)
    }
    
    check_width <- length(unique(width(input))) != 1
    
    if (check_width) {
        
        input_width <- width(input)
        mean_width <- max(input_width)
        which_below <- which(input_width < mean_width)
        howmany <- mean_width - input_width[which_below]
        
        if (type == "fastq") {
            

            to_paste_s <- cpp_create_stringvec(howmany, "N")  # cpp source
            seqs[which_below] <- paste0(seqs[which_below], to_paste_s)
            

            to_paste_q <- cpp_create_stringvec(howmany, "}")
            qual[which_below] <- paste0(qual[which_below], to_paste_q)
            
            out <- ShortReadQ(sread = seqs, quality = qual, id = ids)
            attr(out, "editedWidth") <- TRUE
            return(out)
            
        } else if (type == "sequence") {
            
            to_paste_s <- cpp_create_stringvec(howmany, "N")
            input[which_below] <- paste0(input[which_below], to_paste_s)
            
            attr(input, "editedWidth") <- TRUE
            return(input)
            
        } else if (type == "quality") {
            to_paste_q <- cpp_create_stringvec(howmany, "}")
            input[which_below] <- paste0(input[which_below], to_paste_q)
            
            attr(input, "editedWidth") <- TRUE
            return(input)
            
        }

    } else {
        attr(input, "editedWidth") <- FALSE
        return(input)
    }
}


#' Inject a letter in a set of sequences at random positions 
#' @param my_seq character vector with sequences to inject
#' @param how_many_seqs How many sequences pick to inject Ns. 
#' An interval [min_s, max_s] with min_s minimum 
#' and max_s maximum sequences can be passed. In this case, 
#' a value is picked from the interval. 
#' If NULL, a random value within the interval [1, length(my_seq)] is picked.  
#' @param how_many_letters How many times inject the letter 
#' in the i sequences that are going to be injected. An interval [min_i max_i]
#' can be passed. In this case, a value is randomly 
#' picked for each sequence i. This value represents the number of times
#' that the letter will be injected in the sequence i.
#' If NULL, a random value within the interval [1, width(my_seq[i])] 
#' is picked for each sequence i.
#' @param letter Letter to inject. Default: 'N'
#' @param seed Random seed. Default: NULL
#' @return  character vector
#' @examples 
#' 
#' 
#' s <- random_seq(slength = 10, swidth = 20, seed = 10)
#' s <- inject_letter_random(s, how_many_seqs = 1:30, how_many= 2:10, seed = 10)
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export


inject_letter_random <- function(my_seq, how_many_seqs = NULL, 
                                how_many_letters = NULL, 
    seed = NULL, letter = "N") {
    
    if (!is.null(seed)) 
        set.seed(seed)
    
    if (!is.null(how_many_seqs)) {
        if (any(how_many_seqs < 0)) 
            stop("how_many_seqs must be > 0")
    }
    
    if (!is.null(how_many_letters)) {
        if (any(how_many_letters < 0)) 
            stop("how_many_letters must be > 0")
    }
    
    
    if (is.null(how_many_seqs)) {
        how_many_seqs <- sample(length(my_seq), 1)
    } else {
        how_many_seqs <- c(max(min(how_many_seqs), 1), 
        min(max(how_many_seqs), length(my_seq)))
        how_many_seqs <- sample(how_many_seqs[1]:how_many_seqs[2], 1)
    }
    
    inject_to <- sample(seq_along(my_seq), size = how_many_seqs)
    
    for (i in inject_to) {
        this_seq <- strsplit(as.character(my_seq[i]), "")[[1]]
        this_seq_len <- length(this_seq)
        
        if (is.null(how_many_letters)) {
            this_how_many <- sample(seq_len(this_seq_len), 1)
        } else {
            this_how_many <- c(max(min(how_many_letters), 1), 
            min(max(how_many_letters), this_seq_len))
            this_how_many <- sample(this_how_many[1]:this_how_many[2], 1)
        }
        
        positions <- sample(this_seq_len, size = this_how_many)
        this_seq[positions] <- letter
        this_seq <- paste(this_seq, collapse = "")
        my_seq[i] <- this_seq
    }
    my_seq
}


#' Create sequences names
#' @param n Number of reads
#' @param base_name Base name for strings
#' @param sep Character separing base names and the read number. Default: '_
#' @description Create \code{\link[Biostrings:XStringSet-class]{BStringSet}} 
#' object with names
#' @return  \code{\link[Biostrings:XStringSet-class]{BStringSet}} object
#' @examples 
#' snames <- seq_names(10)
#' snames
#' snames2 <- seq_names(10, base_name = 's', sep = '.')
#' snames2
#' @export

seq_names <- function(n, base_name = "s", sep = "_") {
    BStringSet(paste0(base_name, sep, seq_len(n)))
}

#' Create a named object with random sequences and qualities
#' @param n number of sequences
#' @param widths width of the sequences
#' @param random_widths width must be picked at random from the 
#' passed parameter 'widths', considering the value
#' as an interval where any integer can be picked. Default TRUE. 
#' Otherwise, widths are picked only from the vector passed.
#' @param replace sample widths with replacement? Default TRUE.
#' @param nuc create sequences of DNA (nucleotides = c('A', 'C', 'G', 'T')) or
#'  RNA (nucleotides = c('A, 'C', 'G', 'U'))?. Default: 'DNA'
#' @param len_prob vector with probabilities for each width value. 
#' Default NULL (equiprobability)
#' @param seq_prob a vector of four probabilities values to set 
#' the frequency of the nucleotides 'A', 'C', 'G', 'T', for DNA,
#' or 'A', 'C', 'G', 'U', for RNA. For example = c(0.25, 0.25, 0.5, 0).
#' Default is = c(0.25, 0.25, 0.25, 0.25) (equiprobability for the 4 bases).
#' If the sum of the probabilities is > 1, the values will be nomalized 
#' to the range [0, 1].
#' @param q_prob a vector of range = range(qual), with probabilities to set 
#' the frequency of each quality value. Default is equiprobability.
#' If the sum of the probabilities is > 1, the values will be nomalized 
#' to the range [0, 1].
#' @param qual quality range for the sequences. It must be a range included 
#' in the selected encoding:
#' 
#' 'Sanger' = [0, 40]
#' 
#' 'Illumina1.8' = [0, 41]
#' 
#' 'Illumina1.5' = [0, 40]
#' 
#' 'Illumina1.3' = [3, 40]
#' 
#' 'Solexa' = [-5, 40]
#'
#' example: for a range from 20 to 30 in Sanger encoding, 
#' pass the argument = c(20, 30)
#' @param encod sequence encoding
#' @param base_name Base name for strings
#' @param sep Character separing base names and the read number. Default: '_'
#' @param seed random seed. Default: NULL
#' @return  \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}} object
#' @description Create a \code{\link[ShortRead:ShortReadQ-class]{ShortReadQ}}
#' object with random sequences and qualities
#' @examples 
#' 
#' s1 <- random_seq(slength = 10, swidth = 20, seed = 10)
#' s1
#' 
#' s2 <- random_seq(slength = 10, swidth = 20, seed = 10, 
#' prob = c(0.6, 0.1, 0.3, 0))
#' s2
#' 
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @export

random_length <- function(n, widths, random_widths = TRUE,
                        replace = TRUE, len_prob = NULL, 
                        seq_prob = c(0.25, 0.25, 0.25, 0.25),
                        q_prob = NULL, nuc = c("DNA", "RNA"), qual = NULL, 
                        encod = c("Sanger", "Illumina1.8", "Illumina1.5",
                                "Illumina1.3", "Solexa"), base_name = "s", 
                        sep = "_", seed = NULL) {
    
    if (!is.null(seed)) 
        set.seed(seed)
    
    if (any(widths < 1)) {
        stop("lengths must be > 1")
    }
    if (random_widths) {
        widths <- sample(min(widths):max(widths), n, replace = replace)
    } else {
        widths <- sample(widths, n, replace = replace, prob = len_prob)
    }
    
    s_out <- random_seq(slength = n, swidth = widths, nuc = nuc, 
                        seed = seed, prob = seq_prob)
    
    
    q_out <- random_qual(slength = n, swidth = widths, 
                        qual = qual, encod = encod, 
                        seed = seed, prob = q_prob)
    
    n_out <- seq_names(n = n, base_name = base_name, sep = sep)
    
    ShortReadQ(s_out, q_out, n_out)
}
