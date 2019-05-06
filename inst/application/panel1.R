
panel1 <- function(input, output, session, messages, rList) {
    
    output$N <- renderPrint({
        cat(messages()$nFilterCatch, sep = "\n")
    })
    
    output$filepath <- renderPrint({
        
        if (any(rList()$filepath == "")) {
            cat("No file selected")
        } else {
            # if(dataInput() == '[1] \'error\'\n') { dataInput <- NULL return('Invalid
            # file') } else {
            cat(rList()$filepath, sep = "\n")
        }
    })
    
    output$lowComplexFilter <- renderPrint({
        cat(messages()$lowComplexFilterCatch, sep = "\n")
    })
    
    output$adaptFilter <- renderPrint({
        cat(messages()$seqFilterCatch, sep = "\n")
    })
    
    output$meanQFilter <- renderPrint({
        if (!is.null(rList()$meanQ) && is.null(rList()$thresmeanQ)) {
            cat("Please select a quality threshold value")
        } else {
            cat(messages()$meanQFilterCatch, sep = "\n")
        }
    })
    
    output$TrimQual <- renderPrint({
        if (!is.null(rList()$rm.qual) && is.null(rList()$thresQual)) {
            cat("Please select a quality threshold value")
        } else {
            cat(messages()$trim3FilterCatch, sep = "\n")
        }
    })
    
    output$trimTails <- renderPrint({
        if (!is.null(rList()$rmFixed) && is.null(rList()$rm3) && is.null(rList()$rm5)) {
            cat("Please select remove from 3' and/or 5'")
        } else {
            cat(messages()$fixedFilterCatch, sep = "\n")
        }
    })
    
    output$sizeFilter <- renderPrint({
        if (!is.null(rList()$rmSize) && is.null(rList()$rmMin) && is.null(rList()$rmMax)) {
            cat("Please select min and/or max values")
        } else {
            cat(messages()$lengthFilterCatch, sep = "\n")
        }
    })
    
    output$DuplicFilter <- renderPrint({
        cat(messages()$duplicFilterCatch, sep = "\n")
    })  # // END SALIDAS VENTANA 2 -->
    
    
    output$outResult <- renderPrint({
        if (any(rList()$filepath == "")) {
            cat("No file processed")
        } else {
            cat(messages()$outResult, sep = "\n")
        }
    })
}
