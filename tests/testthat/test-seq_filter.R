context("test-seq_filter.R")
skip_on_cran()

require("Biostrings")
require("ShortRead")

set.seed(10)
input <- random_length(30, 3:7)

test_that("seq filter works", {
  
rm.seq  = c("GAGA", "ATTTAT", "TAGC")
match_before <- unlist(lapply(rm.seq, function(x) grep(x, as.character(sread(input)))))
filtered <- seq_filter(input,rm.seq =  rm.seq)
match_after <- unlist(lapply(rm.seq, function(x) grep(x, as.character(sread(filtered)))))

expect_gt(length(match_before),  0)
expect_equal(length(match_after), 0)
})
