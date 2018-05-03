context("test-trim3q_filter.R")
skip_on_cran()

require("Biostrings")
require("ShortRead")

# create 6 sequences of width 20

input <- random_seq(6, 20, seed = 10)

# create qualities of width 15 and paste to qualities of length 5 used for the tails.
# for two of the sequences, put low qualities in tails

input_q <- random_qual(c(30,40), slength = 6, swidth = 15, seed = 10, encod = "Sanger")
tails <-   random_qual(c(30,40), slength = 6, swidth = 5, seed = 10, encod = "Sanger")
tails[2:3] <- random_qual(c(3, 20), slength = 2, swidth = 5, seed = 10, encod = "Sanger")
input_q <- paste0(input_q, tails)
input_q <- BStringSet(input_q)

# create names
input_names <- seq_names(6)

# create ShortReadQ object
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)


test_that("trim3q filter works", {
  
  # apply the filter 
  filtered <- trim3q_filter(my_read, rm.3qual = 28)

  expect_true(all(width(filtered)[-c(2,3)] == 20))
  expect_true(all(width(filtered)[c(2,3)] == 15))
})

