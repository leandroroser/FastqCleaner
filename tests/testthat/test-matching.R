context("test-matching.R")
skip_on_cran()

library(Biostrings)
Lpattern <- Rpattern <- "TCATG"
subject <- DNAStringSet(c("ATCATGCCATCATGAT", "ATTTGGAATCAT",
"CATGATATTA", "TCATG", "AAAAAA", "AGGTCATG"))

test_that("cutLseq works", {

  # cutLseq
  out <- cutLseq(subject, Lpattern, anchored = FALSE, ranges = FALSE,  min_match_flank = 1)
  expect_equal(width(out), c(2L, 12L, 6L, 0L, 6L, 0L))

  out <- cutLseq(subject, Lpattern, anchored = TRUE, ranges = FALSE, min_match_flank = 1)
  expect_equal(width(out), c(16L, 12L, 6L, 0L, 6L, 8L))
  
  
})


test_that("cutRseq works", {
  
  out <- FastqCleaner:::cutRseq(subject, Rpattern, anchored = FALSE, ranges = FALSE)
  expect_equal(width(out),  c(1L, 8L, 10L, 0L, 6L, 3L))
  
  out <- FastqCleaner:::cutRseq(subject, Rpattern, anchored = TRUE, ranges = FALSE)
  expect_equal(width(out),  c(16L, 8L, 10L, 0L, 6L, 3L))
  
})
