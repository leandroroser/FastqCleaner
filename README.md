<span><h1>  FastqCleaner [![Build Status](https://travis-ci.org/leandroroser/FastqCleaner.svg?branch=master)](https://travis-ci.org/leandroroser/FastqCleaner) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/leandroroser/FastqCleaner?branch=master&svg=true)](https://ci.appveyor.com/project/leandroroser/FastqCleaner) [![R](https://img.shields.io/badge/R%3E%3D-3.0-red.svg)]() [![CRAN](https://img.shields.io/cran/l/devtools.svg)]()  </h1> </span>

A Shiny web-app to clean FASTQ files with R and Bioconductor



## 1. Prerequisites and dependences for the direct installation of this GitHub repository

#### Linux users

- Libraries OpenSSL and libcurl installed

-> Ubuntu:

Paste the following command in terminal:

```diff
sudo apt-get install libssl-dev libcurl4-openssl-dev
```

#### Mac users

- Command Line Tools

Type in a Terminal:
```
xcode-select --install
```
A software update popup will ask if the Command Line Tools must be installed

- XQuartz 

XQuartz can be downloaded from the following link: https://www.xquartz.org/

#### Windows users

- Rtools 

The program is available at https://cran.r-project.org/bin/windows/Rtools/


## 2. Installation

Paste the following command in the R console:
<div style="font-family:courier, courier new, serif;">

```diff
source("https://goo.gl/u2xraS")
```

## 3. Launching the application

```diff
library("FastqCleaner")
launch_fqc()
```

## 4. Description and usage
See the tutorial of the application in <a href = https://leandroroser.github.io/FastqCleaner/> this link </a>
