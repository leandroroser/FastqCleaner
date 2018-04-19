context("test-unique_filter.R")

require("Biostrings")
require("ShortRead")

set.seed(10)
s <- random_seq(10, 10)
s <- sample(s, 30, replace = TRUE)
q <- random_qual(30, 10)
n <- seq_names(30)

my_read <- ShortReadQ(sread = s, quality = q, id = n)

test_that("unique filter works", {
  # apply the filter
filtered <- unique_filter(my_read)
expect_true(!all(isUnique(as.character(sread(my_read)))))  
expect_true(all(isUnique(as.character(sread(filtered)))))  
})

