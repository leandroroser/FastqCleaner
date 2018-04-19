context("test-matching.R")

library(Biostrings)
Lpattern <- Rpattern <- "TCATG"
subject <- DNAStringSet(c("ATCATGCCATCATGAT", "ATTTGGAATCAT",
"CATGATATTA", "TCATG", "AAAAAA", "AGGTCATG"))

test_that("cutLseq works", {

  # cutLseq
  out <- cutLseq(subject, Lpattern, anchored = FALSE, ranges = FALSE)
  expect_equal(width(out), c(2L, 12L, 6L, 0L, 6L, 0L))

  out <- cutLseq(subject, Lpattern, anchored = TRUE, ranges = FALSE)
  expect_equal(width(out), c(16L, 12L, 6L, 0L, 6L, 8L))
  
  # cutLseq(subject, Lpattern, method = "er", error_rate = 0.2, ranges = FALSE)
  # cutLseq(subject, Lpattern, method = "er", error_rate = 0.2,
  #         with.indels = TRUE, ranges = FALSE)
  
})


test_that("cutRseq works", {
  
  out <- cutRseq(subject, Rpattern, anchored = FALSE, ranges = FALSE)
  expect_equal(width(out),  c(1L, 8L, 10L, 0L, 6L, 3L))
  
  out <- cutRseq(subject, Rpattern, anchored = TRUE, ranges = FALSE)
  expect_equal(width(out),  c(15L, 8L, 10L, 0L, 6L, 3L))
  
})
