<span><h1>  FastqCleaner [![Build Status](https://travis-ci.org/leandroroser/FastqCleaner.svg?branch=master)](https://travis-ci.org/leandroroser/FastqCleaner) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/leandroroser/FastqCleaner?branch=master&svg=true)](https://ci.appveyor.com/project/leandroroser/FastqCleaner) [![R](https://img.shields.io/badge/R%3E%3D-3.0-red.svg)]() [![CRAN](https://img.shields.io/cran/l/devtools.svg)]()  </h1> </span>

A Shiny web-app to clean FASTQ files with R and Bioconductor



## This is a development repository
### Please install the application using the canonical way using Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("FastqCleaner")


## Note: due to a change in Shiny, the Bioconductor version is working with Shiny <= 1.6. For the development version, please install using the following code:

install.packages("devtools")
devtools::install_github("https://github.com/leandroroser/FastqCleaner")

- The development version will be pushed to Bioconductor after proper checks (note from the author, 2022-07-01)


