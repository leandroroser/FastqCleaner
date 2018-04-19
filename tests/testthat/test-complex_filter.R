context("test-complex_filter.R")

require("Biostrings")
require("ShortRead")

# create  sequences of different width

input <- lapply(c(0, 6, 10, 16, 20, 26, 30, 36, 40), 
                function(x) random_seq(1, x, seed = 10))


# create repetitive "CG" sequences with length adequante 
# for a total length:
# input +  CG = 40

CG <- lapply(c(20, 17, 15, 12, 10, 7, 5, 2, 0), 
             function(x) paste(rep("CG", x), collapse = ""))


# concatenate input and CG
input  <- mapply("paste", input, CG, sep = "")
input <- DNAStringSet(input)


# create qualities of widths 40

input_q <- random_qual(c(30,40), slength = 9, swidth = 40, 
                       seed = 10, encod = "Sanger")


# create names
input_names <- seq_names(9)

# create ShortReadQ object
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)

test_that("complex filter works", {

  # apply the filter,
  filtered <- complex_filter(my_read)
  
  freq <- dinucleotideFrequency(sread(filtered))
  freq  <- freq /rowSums(freq)
  H <- -rowSums(freq  * log2(freq), na.rm = TRUE)
  H_max <- 3.908135  # max entropy
  
  expect_that(max(H), is_less_than(H_max))

})
