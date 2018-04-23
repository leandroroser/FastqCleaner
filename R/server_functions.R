
#' check onclick
#' @description Function to put a tickmark on click
#' @return Change value of reactive output, without return 
#' @keywords internal

check_onclick_ <- function(.menu_react, .butt_number, my_envir) {

output <- NULL

check_onclick <- function(menu_react, butt_number) {
element <- paste0("check", butt_number)
output[[element]] <- shiny::renderPrint({
    if (menu_react() == TRUE) {
        shiny::HTML(paste0("<script>
        \t\t\t\t\t   $(\"#button", 
        butt_number, " span\").addClass(\"checked\");
        \t\t\t\t\t   </script>"))
    } else {
        shiny::HTML(paste0("<script>
        \t\t\t\t\t   if($(\"#button", 
        butt_number, " span\").hasClass(\"checked\")) {
        \t\t\t\t\t   $(\"#button", 
        butt_number, " span\").removeClass(\"checked\");}
        \t\t\t\t\t   
        \t\t\t\t\t   </script>"))
    }
})
}
environment(check_onclick) <- my_envir
check_onclick(.menu_react, .butt_number)
}

#' outputClean_
#' @return Vector with chunks length and width  information 
#' @keywords internal

outputClean_ <- function(.myFile, .lengthWidthVec, my_envir) {

outputClean <- function(myFile, lengthWidthVec) {
file_length <- length(myFile)
file_width <- Biostrings::width(ShortRead::sread(myFile))
lengthWidthVec[[1]] <- sum(c(lengthWidthVec[[1]], file_length))
lengthWidthVec[[2]] <- unique(c(lengthWidthVec[[2]], file_width))
if (length(lengthWidthVec[[2]]) > 1) {
    lengthWidthVec[[2]] <- c(min(lengthWidthVec[[2]]), 
    max(lengthWidthVec[[2]]))
}
lengthWidthVec
}
environment(outputClean) <- my_envir
outputClean(.myFile, .lengthWidthVec)
}

#' messageFun_
#' @return Changes the state of reactive vector, without return 
#' @keywords internal


messageFun_ <- function(.who, .chunck, .which_read, my_envir) {
messageFun <- function(who, chunck, which_read, envir) {

readsWidth[[who]][[which_read]] <- outputClean_(chunck, 
readsWidth[[who]][[which_read]], envir)
chunck_length <- readsWidth[[who]][[which_read]][[1]]
chunck_width <- readsWidth[[who]][[which_read]][[2]]
if (length(chunck_width) > 1) {
    chunck_width <- paste0(min(chunck_width), "-", 
    max(chunck_width), " cicles")
} else {
    chunck_width <- paste0(readsWidth[[who]][[which_read]][[2]], 
    " cicles")
}
messages[[who]][which_read] <- paste0("Processed: ", chunck_length, 
                                        " reads. ", 
                                        chunck_width)
}
environment(messageFun) <- my_envir
messageFun(.who, .chunck, .which_read, my_envir)
}


#' create_cleanfunction_
#' @description Create a function to process FASTQ files in function
#' of the Shiny parameters selected by the user
#' @return Function with selected cleaning operations
#' @keywords internal

create_cleanfunction_ <- function(my_envir, 
.which_read = c("FORWARD", "REVERSE")) {
.which_read <- match.arg(.which_read)
force(my_envir)  

out <- function() {}

NCleanData <- lowComplexData <- SeqInput <- function() {}
FilterbymeanQ <- TrimCleanData <- FixedCleanData <- function() {}
SizeCleanData <- FilterbymeanQ <- TrimCleanData <- function() {}
FixedCleanData <- SizeCleanData <- DuplicCleanData <- function() {}

fun_body <- function(which_read) {

if (NCleanData()) {
    e1 <- expression(fileToProcess <- n_filter(fileToProcess,
    as.numeric(input$rm.N)), 
    messageFun_("nFilterCatch", 
    fileToProcess, read, my_envir))
} else {
    e1 <- expression()
}

if (lowComplexData()) {
    e2 <- expression(fileToProcess <- complex_filter(fileToProcess,
    input$complexThres), 
    messageFun_("lowComplexFilterCatch", 
    fileToProcess, read, my_envir))
} else {
    e2 <- expression()
}

if (SeqInput()) {
    if (which_read == "FORWARD") {
        e3 <- expression(fileToProcess <- adapter_filter(fileToProcess,
            Lpattern = input$LpatternF, 
            Rpattern = input$RpatternF,
            method = input$adapt_method,
            rc.L = input$reverseLF, 
            rc.R = input$reverseRF, 
            first = input$first_forward,
            with_indels = input$indels, 
            error_rate = input$errate, 
            anchored = input$anchored, 
            min_match_flank = input$min_match_flank, 
            checks = FALSE), 
        messageFun_("seqFilterCatch", fileToProcess, 
        .which_read = 1,  my_envir))
} else {
        e3 <- expression(fileToProcess <- adapter_filter(fileToProcess, 
            Lpattern = input$LpatternR, 
            Rpattern = input$RpatternR, 
            method = input$adapt_method,
            rc.L = input$reverseLR, 
            rc.R = input$reverseRR,
            first = input$first_reverse, 
            with_indels = input$indels, 
            error_rate = input$errate,
            anchored = input$anchored,
            min_match_flank = input$min_match_flank, 
            checks = FALSE), 
        messageFun_("seqFilterCatch", fileToProcess, 
        .which_read = 2,  my_envir))
}  # end 'REVERSE'

} else {
    e3 <- expression()
}

if (FilterbymeanQ()) {
    e4 <- expression(fileToProcess <- qmean_filter(fileToProcess,
        input$thresmeanQ, 
        q_format = my_encoding[["value"]]), 
    messageFun_("meanQFilterCatch", fileToProcess, read, my_envir))
} else {
e4 <- expression()
}

if (TrimCleanData()) {
    e5 <- expression(fileToProcess <- trim3q_filter(fileToProcess, 
        input$thresQual, 
        q_format = my_encoding[["value"]]), 
        messageFun_("trim3FilterCatch", fileToProcess, read, my_envir))
} else {
e5 <- expression()
}

if (FixedCleanData()) {
    e6 <- expression(fileToProcess <- fixed_filter(fileToProcess, input$rm3, 
    input$rm5), 
    messageFun_("fixedFilterCatch", fileToProcess, read, my_envir))
} else {
    e6 <- expression()
}

if (SizeCleanData()) {
    e7 <- expression(fileToProcess <- length_filter(fileToProcess, input$rmMin, 
    input$rmMax), 
    messageFun_("lengthFilterCatch", fileToProcess, read,  my_envir))
} else {
e7 <- expression()
}

if (DuplicCleanData()) {
    e8 <- expression(fileToProcess <- unique_filter(fileToProcess),
    messageFun_("duplicFilterCatch", fileToProcess, read, my_envir))
} else {
    e8 <- expression()
}
    e9 <- expression(fileToProcess)
    as.call(c(as.name("{"), e1, e2, e3, e4, e5, e6, e7, e8, e9))
}

environment(fun_body) <- my_envir
body(out) <- fun_body(which_read = .which_read)
formals(out) <- alist(fileToProcess = , read = , my_envir = )
environment(out) <- my_envir

out
}

#' processingFunction_
#' @description This function is the core of the application. It is used 
#' for the program to process the FASTQ file/s in the environment 
#' of the Shiny app. Note that this program makes a call to
#' create_cleanfunction
#' @return Processes the input FASTQ file, without return
#' @keywords internal


processingFunction_ <- function(my_envir) {

input <- this_envir <- filepath <- session <- NULL

create_processing_function <- function() {

if (input$fileTypeOut == "fastq") {
    compress <- FALSE
} else {
    compress <- TRUE
}


i <- 0


if (input$fileTypeIn == "SR") {
    .cleanfunction <- create_cleanfunction_(this_envir, 
    .which_read = "FORWARD")
    no_output <- as.character(body(.cleanfunction))[2] == "fileToProcess"

} else {
    .cleanfunctionF <- create_cleanfunction_(this_envir, 
    .which_read = "FORWARD")
    .cleanfunctionR <- create_cleanfunction_(this_envir,
    .which_read = "REVERSE")
    no_output <- as.character(body(.cleanfunctionF))[2] == "fileToProcess"


}


if (input$fileTypeIn == "SR" && no_output) {
    outname <- sprintf("%s_out", filepath$x)
if (file.exists(outname)) {
    file.remove(outname)
}
    file.symlink(filepath$x, outname)
}

if (input$fileTypeIn == "PE" && no_output) {

    outname1 <- sprintf("%s_out", (filepath$x)[1])
    outname2 <- sprintf("%s_out", (filepath$x)[2])

if (file.exists(outname1)) {
    file.remove(outname1)
}

if (file.exists(outname2)) {
    file.remove(outname2)
}

file.symlink((filepath$x)[1], outname1)
file.symlink((filepath$x)[2], outname2)
}


if (input$fileTypeIn == "SR" && !no_output) {


if (file.exists(sprintf("%s_out", filepath$x))) {
    file.remove(sprintf("%s_out", filepath$x))
    file.create(sprintf("%s_out", filepath$x))
}

cat("Counting lines in file ...\n")
lfile <- ShortRead::countLines(filepath$x)/4
cat(lfile, "reads found\n\n")
largo <- 0L
ancho <- c(-1L, -1L)


progress <- shiny::Progress$new(session, min = 0, max = lfile)
progress$set(message = "Processing the file",
detail = "This may take a while...")

stream <- ShortRead::FastqStreamer(filepath$x, n = input$nfile)


while (TRUE) {
    passfile <- try(ShortRead::yield(stream), silent = TRUE)
    if(length(passfile) == 0) {

    if(length(unique(ancho)) > 1) {
        cycle_info <- paste0(ancho[1], "-", ancho[2], " cycles")
    } else {
        cycle_info <- paste0(unique(ancho), " cycles")
    }

messages[["outResult"]] <- paste0(largo, " reads. ", cycle_info, 
"\n",
sprintf("%s_out", (filepath$x)))

break

} else {
progress$set(value = i)
cat(i, " reads processed", "\n")
i <- i + length(passfile)

outFile <- .cleanfunction(passfile, read = 1, this_envir)
ShortRead::writeFastq(outFile, sprintf("%s_out", 
filepath$x), mode = "a", 
compress = compress)

largo <- largo + length(outFile)
ancho <- unique(range(c(ancho[ancho != -1], Biostrings::width(outFile))))
cat(ancho, "\n")

}
}

close(stream)
progress$close()
}

if (input$fileTypeIn == "PE" && !no_output) {

    cat("Counting lines in file ...\n")
    lfile <- (ShortRead::countLines((filepath$x)[1]))/4
    cat(lfile, "reads found\n\n")

    largoL <- largoR <- 0L
    anchoL <- anchoR <- c(-1L, -1L)

    for (j in filepath$x) {
        if (file.exists(sprintf("%s_out", j))) {
        file.remove(sprintf("%s_out", j))
    }
        file.create(sprintf("%s_out", j))
    }

    progress <- shiny::Progress$new(session, min = 0, max = lfile)
    progress$set(message = "Processing the file",
    detail = "This may take a while...")

    streamL <- ShortRead::FastqStreamer((filepath$x)[1], n = input$nfile)
    streamR <- ShortRead::FastqStreamer((filepath$x)[2], n = input$nfile)

    while (TRUE) {

        passfileL <- try(ShortRead::yield(streamL), silent = TRUE)
        passfileR <- try(ShortRead::yield(streamR), silent = TRUE)

        if(length(passfileL) == 0 && length(passfileR) == 0) {

            if(length(unique(anchoL)) > 1) {
                cycle_info_L <- paste0(anchoL[1], "-", anchoL[2], " cycles")
            } else {
                cycle_info_L <- paste0(unique(anchoL), " cycles")
            }

            messages[["outResult"]][1] <- paste0(largoL, " reads. ", 
                                            cycle_info_L, "\n",
            sprintf("%s_out", (filepath$x)[1]))

            if (length(unique(anchoR)) > 1) {
                cycle_info_R <- paste0(anchoR[1], "-", anchoR[2], " cycles")
            } else {
                cycle_info_R <- paste0(unique(anchoR), " cycles")
            }

            messages[["outResult"]][2] <- paste0(largoR, " reads. ", 
                                                cycle_info_R, "\n", 
            sprintf("%s_out", (filepath$x)[2]))

            break

        } else {
        progress$set(value = i)
        cat(i, " reads processed", "\n")
        i <- i + length(passfileL)

        outFileL <- .cleanfunctionF(passfileL, read = 1, this_envir)
        outFileR <- .cleanfunctionR(passfileR, read = 2, this_envir)

        inL <- ShortRead::id(passfileL)
        outL <- ShortRead::id(outFileL)
        indexL <- as.character(inL) %in% as.character(outL)
        names(indexL) <- seq(along = indexL)
        indexL <- names(indexL[indexL])

        inR <- ShortRead::id(passfileR)
        outR <- ShortRead::id(outFileR)
        indexL <- as.character(inL) %in% as.character(outL)
        names(indexR) <- seq(along = indexR)
        indexR <- names(indexR[indexR])

        indexpairL <- indexL %in% indexR
        indexpairR <- indexR %in% indexL

        outFileL <- outFileL[indexpairL]
        outFileR <- outFileR[indexpairR]

        ShortRead::writeFastq(outFileL, sprintf("%s_out", (filepath$x)[1]), 
        mode = "a", compress = compress)

        ShortRead::writeFastq(outFileR, sprintf("%s_out", (filepath$x)[2]), 
        mode = "a", compress = compress)


        largoL <- largoL + length(outFileL)
        largoR <- largoR + length(outFileR)

        anchoL <- unique(range(c(anchoL[anchoL != -1], 
                Biostrings::width(outFileL))))
        anchoR <- unique(range(c(anchoR[anchoR != -1], 
                Biostrings::width(outFileR))))

    }
    }

        close(streamL)
        close(streamR)
        progress$close()

        }
    }

environment(create_processing_function) <- my_envir
create_processing_function()
}
