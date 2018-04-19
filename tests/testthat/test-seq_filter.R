context("test-seq_filter.R")

require("Biostrings")
require("ShortRead")

input <- random_length(30, 3:7, seed = 10)

test_that("seq filter works", {
  
rm.seq  = c("TGGTC", "CGGT", "GTTCT", "ATA")
match_before <- unlist(lapply(rm.seq, function(x) grep(x, as.character(sread(input)))))
filtered <- seq_filter(input,rm.seq =  rm.seq)
match_after <- unlist(lapply(rm.seq, function(x) grep(x, as.character(sread(filtered)))))

expect_gt(length(match_before),  0)
expect_equal(length(match_after), 0)
})
