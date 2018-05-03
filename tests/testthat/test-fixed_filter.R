context("test-fixed_filter.R")
skip_on_cran()

require("Biostrings")
require("ShortRead")

# create 6 sequences of width 20

input <- random_seq(6, 20, seed = 10)

# create qualities of width 20 

input_q <- random_qual(c(30,40), slength = 6, swidth = 20, 
                       seed = 10, encod = "Sanger")


# create names
input_names <- seq_names(6)

# create ShortReadQ object
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)


test_that("fixed filter works", {
  # apply the filter 
  filtered3 <- fixed_filter(my_read, trim5 = 5)
  
  filtered5 <- fixed_filter(my_read, trim3 = 5)
  
  filtered3and5 <- fixed_filter(my_read, trim3 = 10, trim5 = 5)
  
  expect_that(unique(width(filtered3)), equals(15))
  expect_that(unique(width(filtered5)), equals(15))
  expect_that(unique(width(filtered3and5)), equals(5))
  
})
