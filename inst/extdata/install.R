
##########################
# FastqCleaner installer
##########################

message("\n")
message("######################")
message("   FastqCleaner")
message("######################\n")
message("Installing dependencies, this may take a while...")
message("\n")

R_repo <- "https://cloud.r-project.org/"

# set internet options
options(download.file.method = "libcurl", url.method = "libcurl")

# set library -------------------------------------------------------------------------------------
where <- normalizePath(Sys.getenv("R_LIBS_USER"))
if(!dir.exists(where)) {
  dir.create(where, recursive =  TRUE)
}
.libPaths(c(.libPaths(), where))

github_where <- paste("--library=", paste0("'", where, "'"))

# Check OS ---------------------------------------------------------------------------------------
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

this_os <- get_os()

# Check dependences for each os

if (this_os == "osx") {
  test_xcode <- class(try(system2("gcc --version", stdout = TRUE), silent = TRUE))
  if(is(test_xcode, "try-error")) {
    
    stop(toupper("Command Line Tools is required. Please paste install Command Line Tools pasting the following 
          command in a Terminal and then selecting install in the popup dialog:\n 
          xcode-select --install"))
  }
  if(length(list.files("/opt/X11/bin", pattern = "Xquartz")) == 0) {
  stop(toupper("X11 is required. Please visit www.xquartz.org to download and install Xquartz."))
  }
}

if(this_os == "linux") {
  is_ubuntu <- grep("ubuntu", sessionInfo()$running,  ignore.case = TRUE)
  if(length(is_ubuntu) > 0) {
    p1 <- is(try(system2("dpkg -l |grep libssl-dev", stdout =  TRUE)), 
                "try-error")
    p2 <- is(try(system2("dpkg -l |grep libcurl4-openssl-dev", stdout =  TRUE)),
                "try-error")
    msg <- ""
    what <- ""
    if(p1) {
      msg <- ("libssl-dev in required\n")
      what <- "libssl-dev "
    }
    if(p2) {
      msg <- c(msg, "libcurl4-openssl-dev in required\n")
      what <- paste0(what, "libcurl4-openssl-dev ")
    }
    if(p1 || p2) {
      stop(toupper(paste0(msg, "Please install the dependences with the following command\n: sudo apt get install ", what)))
    }
  } else {
    cat(toupper("-----------------------------------------------------------------------------------\n
Important: please check that the 'libcurl' and 'libssl' packages corresponding 
to your Linux distribution are installed. If the installation fails, 
                    install them and try again.\n
-----------------------------------------------------------------------------------\n\n"))
  }
  
}

if(this_os == "windows") {
  if(!require(devtools)) {
  cat(toupper("Installing devtools...\n"))
  try(install.packages("devtools", lib = where, repos = R_repo, quiet = TRUE, verbose = FALSE), silent = TRUE)
  }
  
  if(!require(devtools)) {
    warning(toupper("Cannot install devtools. Trying to install Rtools..."))
    install.packages("installr",lib = where, repos = R_repo)
    rtoolsexit <- installr::install.Rtools()
    if(!rtoolsexit) {
      stop("Failed.\nRtools needs to be manually installed from: https://cran.r-project.org/bin/windows/Rtools/")
    }
  }
}


# install dependences------------------------------------------------------------------------------
packages <- c("devtools", "Rcpp", "DT", "shiny", "methods", "graphics", "shinyBS", "RCurl")

for(i in packages) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, lib = where, repos = R_repo)
  }
}

# install Bioconductor dependences----------------------------------------------------------------
bpackages <- c("IRanges", "Biostrings", "ShortRead", "S4Vectors")
<<<<<<< HEAD
source("https://bioconductor.org/biocLite.R")

for(i in bpackages) {
  if (!require(i, character.only = TRUE)) {
    biocLite(i, lib = where,suppressUpdates = TRUE, dependencies = TRUE)
=======
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

for(i in bpackages) {
  if (!require(i, character.only = TRUE)) {
    BiocManager::install(i, lib = where,suppressUpdates = TRUE, dependencies = TRUE)
>>>>>>> upstream/master
  }
}

# install package --------------------------------------------------------------------------------
library("devtools")
withr::with_libpaths(where, 
install_github("leandroroser/FastqCleaner", build_vignettes = TRUE, verbose = TRUE, args = github_where)
)
# EOF

