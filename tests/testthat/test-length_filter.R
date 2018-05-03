context("test-length_filter.R")
skip_on_cran()

require("Biostrings")
require("ShortRead")

#create  ShortReadQ object width widths between 1 and 100
  
input <- random_length(100, widths = seq_len(100), seed = 10)

test_that("length filter works", {
  
  ## apply the filter, removing sequences with  10>length> 80
  filtered <- length_filter(input, rm.min = 10, rm.max = 80)
  
  expect_true(all(width(filtered) >= 10))
  expect_true(all(width(filtered) <= 80))
})
