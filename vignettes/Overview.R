## ----eval=FALSE, include=TRUE----------------------------------------------
#  library('FastqCleaner')
#  launch_fqc()

## ----eval=FALSE, include=TRUE----------------------------------------------
#  download.file("https://goo.gl/hb4Kr9", "example_fastq.gz")

## ----eval=TRUE, echo=TRUE--------------------------------------------------
### Examples

## ----include=FALSE---------------------------------------------------------
require("Biostrings")
require("ShortRead")
require("FastqCleaner")
set.seed(10)

## ----eval = FALSE----------------------------------------------------------
#  require("Biostrings")
#  require("ShortRead")
#  require("FastqCleaner")

## ----eval=TRUE, echo=TRUE--------------------------------------------------

# create sequences
set.seed(10) 
# nota that the use of set.seed before the call to the 
# random generators allows reproducibility of the
# examples

input <- random_seq(6, 43)
input

# create qualities of width 50
set.seed(10)
input_q <- random_qual(c(30,40), slength = 6, swidth = 50, 
encod = "Sanger")

# create names
input_names <- seq_names(length(input))


### FULL ADAPTER IN 3'
adapter <- "ATCGACT"

# Create sequences with adapter
my_seqs <- paste0(input, adapter)
my_seqs <- DNAStringSet(my_seqs)
my_seqs

# create ShortReadQ object
my_read <- ShortReadQ(sread = my_seqs, quality = input_q, id = input_names)

# trim adapter
filtered <- adapter_filter(my_read, Lpattern = adapter)
sread(filtered)

### PARTIAL ADAPTER IN 5'
adapter <- "ATCGACT"
subadapter <- subseq(adapter, 1, 4)

# Create sequences with adapter
my_seqs <- paste0(input, subadapter)
my_seqs <- DNAStringSet(my_seqs)
my_seqs

# create ShortReadQ object
my_read <- ShortReadQ(sread = my_seqs, quality = subseq(input_q, 1, 47), 
id = input_names)

# trim adapter
filtered <- adapter_filter(my_read, Rpattern = adapter)
sread(filtered)


## --------------------------------------------------------------------------

# create  sequences of different width
set.seed(10)
input <- lapply(c(0, 6, 10, 16, 20, 26, 30, 36, 40), 
            function(x) random_seq(1, x))


# create repetitive "CG" sequences with length adequante 
# for a total length input +  CG = 40

CG <- lapply(c(20, 17, 15, 12, 10, 7, 5, 2, 0), 
            function(x) paste(rep("CG", x), collapse = ""))

# concatenate input and CG
input  <- mapply("paste", input, CG, sep = "")
input <- DNAStringSet(input)
input

# plot relative entropy (E, Shannon 1948)
H_plot <- function(x, H_max = 3.908135) {
    freq <- dinucleotideFrequency(x)
    freq  <- freq /rowSums(freq)
    H <- -rowSums(freq  * log2(freq), na.rm = TRUE)
    plot(H/H_max, type="l", xlab = "Sequence", ylab= "E")
    points(H/H_max, col = "#1a81c2", pch = 16, cex = 2)
}

H_plot(input)

## --------------------------------------------------------------------------
# create qualities of widths 40
set.seed(10)
input_q <- random_qual(c(30,40), slength = 9, swidth = 40, 
    encod = "Sanger")

# create names
input_names <- seq_names(9)


# create ShortReadQ object
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)

# apply the filter, 
filtered <- complex_filter(my_read)
sread(filtered)

H_plot(sread(filtered))

## --------------------------------------------------------------------------

# create sequences, qualities and names of width 20
set.seed(10)
input <- random_seq(6, 20)
input

set.seed(10)
input_q <- random_qual(c(30,40), slength = 6, swidth = 20, 
    encod = "Sanger")

input_names <- seq_names(6)

# create ShortReadQ object
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)

# apply the filter 
filtered3 <- fixed_filter(my_read, trim5 = 5)
sread(filtered3)

filtered5 <- fixed_filter(my_read, trim3 = 5)
sread(filtered5)

filtered3and5 <- fixed_filter(my_read, trim3 = 10, trim5 = 5)
sread(filtered3and5)


## --------------------------------------------------------------------------

# create  ShortReadQ object width widths between 1 and 60
set.seed(10)
input <- random_length(10, widths = 1:60)
sread(input)

# apply the filter, removing sequences with  5>length> 30
filtered <- length_filter(input, rm.min = 5, rm.max = 30)
sread(filtered)

## --------------------------------------------------------------------------
# create 10 sequences of width 20
set.seed(10)
input <- random_seq(10, 20)
input

# inject N's
set.seed(10)
input <- inject_letter_random(input, how_many_seqs = 1:5,
    how_many = 1:10)
input

#'  
hist(letterFrequency(input, "N"), breaks = 0:10, 
    main  = "Ns Frequency", xlab = "# Ns",
    col = "#1a81c2")

## --------------------------------------------------------------------------

# Create qualities, names and ShortReadQ object
set.seed(10)
input_q <- random_qual(10, 20)
input_names <- seq_names(10)
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)

# Apply the filter 
filtered <- n_filter(my_read, rm.N = 3)
sread(filtered)
hist(letterFrequency(sread(filtered), "N"), 
    main = "Ns distribution", xlab = "",
    col = "#1a81c2")

## --------------------------------------------------------------------------

# create 30 sequences of width 20, 15 with low quality and 15 with high quality
set.seed(10)
input <- random_seq(30, 20)

set.seed(10)
my_qual_H <- random_qual(c(30,40), slength = 15, swidth = 20,
    encod = "Sanger")

set.seed(10)
my_qual_L <- random_qual(c(5,30), slength = 15, swidth = 20, 
    encod = "Sanger")
input_q<- c(my_qual_H, my_qual_L)

input_names <- seq_names(30)
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)

# Plot of average qualities
qual_plot <- function(x, cutoff) {
q <- alphabetScore(x) / width(x)
plot(q, type="l", xlab = "Sequence", ylab= "Average quality", ylim = c(0, 40))
points(q, col = "#1a81c2", pch = 16, cex = 2)
lines(seq_along(q), rep(cutoff, length(q)), type="l", col = "red", lty=2)
text(length(q), cutoff+2, cutoff)
}

#' Average qualities before
qual_plot(my_read, cutoff = 30)

## --------------------------------------------------------------------------
# Apply the filter
filtered <- qmean_filter(my_read, minq = 30)

# Average qualities after
qual_plot(filtered, cutoff = 30)


## --------------------------------------------------------------------------

# Generate random sequences
set.seed(10)
input <- random_length(30, 3:7)

# Remove sequences that contain the following patterns:
rm.seq  = c("TGGTC", "CGGT", "GTTCT", "ATA")
match_before <- unlist(lapply(rm.seq, function(x) grep(x, 
as.character(sread(input)))))
match_before

filtered <- seq_filter(input,rm.seq =  rm.seq)

# Verify that matching sequences were removed
match_after <- unlist(lapply(rm.seq, function(x) {
                grep(x, as.character(sread(filtered)))}))
match_after

## --------------------------------------------------------------------------

# Create 6 sequences of width 20
set.seed(10)
input <- random_seq(6, 20)
input

# Create Phred+33 qualities of width 15 and paste to qualities of length 
# 5 used for the tails.
# for three of the sequences, put low qualities in tails
set.seed(10)
my_qual <- random_qual(c(30,40), slength = 6, swidth = 15, 
    encod = "Sanger")
set.seed(10)
tails <-   random_qual(c(30,40), slength = 6, swidth = 5, 
    encod = "Sanger")

# Low quality tails in sequences 2, 3 & 4
set.seed(10)
tails[2:4] <- random_qual(c(3, 20), slength = 3, swidth = 5,
    encod = "Sanger")
my_qual <- paste0(my_qual, tails)
input_q <- BStringSet(my_qual)
input_q

# Watch qualities before filtering
as.matrix(PhredQuality(input_q))

# Create names and ShortReadQ object
input_names <- seq_names(6)
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)

# Apply the filter 
filtered <- trim3q_filter(my_read, rm.3qual = 28)
sread(filtered)


## --------------------------------------------------------------------------

# Create duplicated sequences
s <- random_seq(10, 10)
s <- sample(s, 30, replace = TRUE)

# Create a ShortReadQ object
q <- random_qual(30, 10)
n <- seq_names(30)
my_read <- ShortReadQ(sread = s, quality = q, id = n)

# Check presence of duplicates
isUnique(as.character(sread(my_read)))

# Apply the filter
filtered <- unique_filter(my_read)
isUnique(as.character(sread(filtered)))

