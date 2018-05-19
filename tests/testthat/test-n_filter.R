context("test-n_filter.R")
skip_on_cran()

require("Biostrings")
require("ShortRead")
set.seed(10)

# create 6 sequences of width 20
set.seed(10)
input_seq <- random_seq(50, 20)

set.seed(10)
input_seq <- inject_letter_random(input_seq, 30, 
                                  how_many_letters = seq_len(10))
set.seed(10)
input_qual <- random_qual(50, 20)
input_names <- seq_names(50)
my_read <- ShortReadQ(input_seq, input_qual, input_names)

test_that("n_filter works", {
  
# apply the filter 
filtered <- n_filter(my_read, rm.N = 5)

expect_true(all(letterFrequency(sread(filtered), "N") <= 5))

})

