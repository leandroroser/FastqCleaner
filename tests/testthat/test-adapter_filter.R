context("test-adapter_filter.R")

require("Biostrings")
require("ShortRead")
set.seed(10)

input <- random_seq(6, 43, seed = 10)
# create qualities of width 50
input_q <- random_qual(c(30,40), slength = 6, swidth = 50, 
                       seed = 10, encod = "Sanger")

# create names
input_names <- seq_names(length(input))



test_that("test that adapter remotion from right works", {

  ### FULL ADAPTER
  # add adapter in 3'
  adapter <- "ATCGACT"
  ## Remove reverse complement from right
  my_seqs <- paste0(input, adapter)
  my_seqs <- DNAStringSet(my_seqs)
  
  # create ShortReadQ object
  my_read <- ShortReadQ(sread = my_seqs, quality = input_q, id = input_names)
  
  # trim adapter
  filtered <- adapter_filter(my_read, Rpattern = adapter)
  
  expect_equal(unique(width(filtered)), 43)
  
  
  ### PARTIAL ADAPTER
  
  adapter <- "ATCGACT"
  subadapter <- subseq(adapter, 1, 4)
  ## Remove reverse complement from right
  my_seqs <- paste0(input, subadapter)
  my_seqs <- DNAStringSet(my_seqs)
  
  # create ShortReadQ object
  my_read <- ShortReadQ(sread = my_seqs, quality = subseq(input_q, 1, 47), id = input_names)
  
  # trim adapter
  filtered <- adapter_filter(my_read, Rpattern = adapter)
  
  expect_equal(unique(width(filtered)), 43)
  
})


test_that("test partial remotion from right reverse-complemented works", {

  ## TEST RC RIGHT
  
  # add adapter in 3'
  adapter <- "ATCGACT"
  ## Remove reverse complement from right
  my_seqs <- paste0(input, as.character(Biostrings::reverseComplement(DNAString(adapter))))
  my_seqs <- DNAStringSet(my_seqs)
  
# create ShortReadQ object
  my_read <- ShortReadQ(sread = my_seqs, quality = input_q, id = input_names)
  
  # trim adapter
  filtered <- adapter_filter(my_read, Rpattern = adapter, rc.R = TRUE)
 
  expect_equal(unique(width(filtered)), 43)
  
})

test_that("test that partial remotion from left works", {
  
  
  ## FULL ADAPTER
  # add adapter in 3'
  adapter <- "ATCGACT"
  ## Remove reverse complement from right
  my_seqs <- paste0(adapter, input)
  my_seqs <- DNAStringSet(my_seqs)
  
  # create ShortReadQ object
  my_read <- ShortReadQ(sread = my_seqs, quality = input_q, id = input_names)
  
  # trim adapter
  filtered <- adapter_filter(my_read, Lpattern = adapter)
  
  expect_equal(unique(width(filtered)), 43)
  
  
  ### PARTIAL ADAPTER
  adapter <- "ATCGACT"
  subadapter <- subseq(adapter, 4, 7)
  ## Remove reverse complement from right
  my_seqs <- paste0(subadapter, input)
  my_seqs <- DNAStringSet(my_seqs)
  
  # create ShortReadQ object
  my_read <- ShortReadQ(sread = my_seqs, quality = subseq(input_q, 1, 47), id = input_names)
  
  # trim adapter
  filtered <- adapter_filter(my_read, Lpattern = adapter)
  
  expect_equal(unique(width(filtered)), 43)
  
})


test_that("test that partial remotion from left reverse-complemented works", {
  
  # add adapter in 3'
  adapter <- "ATCGACT"
  ## Remove reverse complement from right
  my_seqs <- paste0( as.character(Biostrings::reverseComplement(DNAString(adapter))), input)
  my_seqs <- DNAStringSet(my_seqs)
  
  # create ShortReadQ object
  my_read <- ShortReadQ(sread = my_seqs, quality = input_q, id = input_names)
  
  # trim adapter
  filtered <- adapter_filter(my_read, Lpattern = adapter, rc.L = TRUE)
  
  expect_equal(unique(width(filtered)), 43)
  
})


test_that("test that remotion of partial adapters within sequences works", {
  
## FROM RIGHT 
# create qualities of width 50
input_q <- random_qual(c(30,40), slength = 6, 
                       swidth = 47, seed = 10, encod = "Sanger")

# create names
input_names <- BStringSet(paste0("v", seq_len(length(input))))
adapter <- "ATCGACT"
## Remove reverse complement from right
my_seqs <- paste0(random_seq(6, 20, seed = 10), 
                             adapter, random_seq(6, 20, seed = 10))
my_seqs <- DNAStringSet(my_seqs)

# create ShortReadQ object
my_read <- ShortReadQ(sread = my_seqs, quality = input_q, id = input_names)

# trim adapter
filtered <- adapter_filter(my_read, Rpattern = adapter)

expect_equal(unique(width(filtered[-1])), 47)

### FROM LEFT
filtered <- adapter_filter(my_read, Lpattern = adapter)

expect_equal(unique(width(filtered[-1])), 47)

})

