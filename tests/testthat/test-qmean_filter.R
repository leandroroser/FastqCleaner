context("test-qmean_filter.R")

require("Biostrings")
require("ShortRead")


# create 6 sequences of width 20
input <- random_seq(30, 20, seed = 10)

# create qualities of width 20 
## high quality (15 sequences)
my_qual <- random_qual(c(30,40), slength = 15, swidth = 20,
                       seed = 10, encod = "Sanger")
## low quality (15 sequences)
my_qual_2 <-   random_qual(c(5,30), slength = 15, swidth = 20, 
                           seed = 10, encod = "Sanger")

# concatenate vectors
input_q<- c(my_qual, my_qual_2)

# create names
input_names <- seq_names(30)

# create ShortReadQ object
my_read <- ShortReadQ(sread = input, quality = input_q, id = input_names)


test_that("q_mean filter works", {

  filtered <- qmean_filter(my_read, minq = 30)
  
  # average qualities
  avq <- alphabetScore(filtered) / width(filtered)
  
  expect_true(all(avq >=30))
})
