
check_onclick_ <- FastqCleaner:::check_onclick_
messageFun_ <- FastqCleaner:::messageFun_
processingFunction_ <- FastqCleaner:::processingFunction_
outputClean_ <- FastqCleaner:::outputClean_
create_cleanfunction_ <- FastqCleaner:::create_cleanfunction_
check_encoding <- FastqCleaner::check_encoding


shinyServer(
  function(input, output, session) {
    
    if (Sys.getenv("SHINY_PORT") == "") 
      options(shiny.maxRequestSize = 10000 * 1024^2)
    
    session$onSessionEnded(function() {
      stopApp()
    })  
    this_envir <- environment()
    
    
    hasRun <- reactiveValues(x = FALSE)
    
    
    output$myPath <- renderPrint(cat("choose a file..."))
    output$encoding <- renderPrint(cat("..."))
    
    output$runDone <- renderPrint({
      HTML(paste0(""))
    })
    
    tPlot <- reactiveValues(x = list())
    
    filepath <- reactiveValues(x = "")
    
    observeEvent(input$fileTypeIn, {
      filepath$x <- ""
      output$myPath <- renderPrint(cat("choose a file..."))
      output$encoding <- renderPrint(cat("..."))
      tPlot$x <- list(HTML("<script></script>"))
    })
    
    
    readsWidth <- reactiveValues(nFilterCatch = list(list(numeric(), numeric()), 
      list(numeric(), numeric())), 
      lowComplexFilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())),
      seqFilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())),
      meanQFilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())), 
      trim3FilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())), 
      fixedFilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())), 
      lengthFilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())), 
      duplicFilterCatch = list(list(numeric(), 
        numeric()), list(numeric(), numeric())))
    
    # setear mensajes at start
    messages <- reactiveValues(nFilterCatch = c("Filter unused", ""), 
      lowComplexFilterCatch = c("Filter unused", ""), 
      seqFilterCatch = c("Filter unused", ""), 
      meanQFilterCatch = c("Filter unused", ""), 
      trim3FilterCatch = c("Filter unused", ""), 
      fixedFilterCatch = c("Filter unused", ""), 
      lengthFilterCatch = c("Filter unused", ""), 
      duplicFilterCatch = c("Filter unused", ""), 
      outResult = c("No file processed yet", ""))
    
    
    observe({
      if (input$resetSelections == 0) {
        return()
      }
      isolate({
        
        passfile <- NULL
        vars <- c("nFilterCatch", "lowComplexFilterCatch", "seqFilterCatch", 
          "meanQFilterCatch", "trim3FilterCatch", "fixedFilterCatch", 
          "lengthFilterCatch", "duplicFilterCatch")
        
        for (i in vars) {
          messages[[i]] <- c("Filter unused", "")
          readsWidth[[i]] <- list(list(numeric(), 
            numeric()), list(numeric(), 
              numeric()))
        }
        
        messages$outResult = c("No file processed yet", "")
        
        inputNames <- c("rm.N", "rm.lowcomplex", "rm.adapt", 
          "rm.qual", "rmFixed", 
          "rmSize", "rm.duplic", "nuc")
        
        for (i in inputNames) {
          updateCheckboxGroupInput(session, i, selected = character(0))
        }
        
        output$runDone <- renderPrint({
          HTML(paste0("<script> 
            $(\"#button1 span\").removeClass(\"checked\");
            $(\"#button2 span\").removeClass(\"checked\");
            $(\"#button3 span\").removeClass(\"checked\");
            $(\"#button4 span\").removeClass(\"checked\");
            $(\"#button5 span\").removeClass(\"checked\");
            $(\"#button6 span\").removeClass(\"checked\");
            </script>"))
        })
        
        showModal(modalDialog(title = "FastqCleaner message", 
          "Selections reseted!", 
          footer = modalButton("Ok")))
        
        })
    })
    
    observeEvent(input$resetAll, {
      
      filepath$x <- ""
      output$myPath <- renderPrint(cat("choose a file..."))
      output$encoding <- renderPrint(cat("..."))
      passfile <- NULL
      
      vars <- c("nFilterCatch", "lowComplexFilterCatch", 
        "seqFilterCatch", "meanQFilterCatch",
        "trim3FilterCatch", "fixedFilterCatch", 
        "lengthFilterCatch", "duplicFilterCatch")
      
      for (i in vars) {
        messages[[i]] <- c("Filter unused", "")
        readsWidth[[i]] <- list(list(numeric(), numeric()), 
          list(numeric(), numeric()))
      }
      
      messages$outResult = c("No file processed yet", "")
      
      inputNames <- c("rm.N", "rm.lowcomplex", "rm.adapt",
        "rm.qual", "rmFixed", 
        "rmSize", "rm.duplic", "nuc")
      
      for (i in inputNames) {
        updateCheckboxGroupInput(session, i, selected = character(0))
      }
      output$runDone <- renderPrint({
        HTML(paste0("<script> 
          $(\"#button1 span\").removeClass(\"checked\");
          $(\"#button2 span\").removeClass(\"checked\");
          $(\"#button3 span\").removeClass(\"checked\");
          $(\"#button4 span\").removeClass(\"checked\");
          $(\"#button5 span\").removeClass(\"checked\");
          $(\"#button6 span\").removeClass(\"checked\");
          </script>"))
      })
      
      output$set_encoding_radio <- renderText({
        out <- list()
        for (i in 1:5) {
          out[[i + 1]] <- HTML({
            paste0("<script>
              $( \"input[name='customEncod'][value ='", 
              i, "']\" ).prop(\"checked\", false);
              </script>")
          })
          }
        out[[1]] <- HTML({
          paste0("<script>
            $( \"input[name='customEncod'][value ='", 
            0, "']\" ).prop(\"checked\", true);
            </script>")
        })
        unlist(out)
        })  
      
      showModal(modalDialog(title = "FastqCleaner message",
        "All reseted!", footer = modalButton("Ok")))
      })
    
    
    error_selection <- reactiveValues(error_maxN = FALSE, 
      error_complexThres = FALSE,
      error_errate = FALSE, 
      error_thresmeanQ = FALSE,
      error_thresQual = FALSE, 
      error_rm35 = FALSE, 
      error_rmMinMax = FALSE)
    
    
    
    NCleanData <- reactive({
      if (!is.null(input$rm.N)) {
        TRUE
      } else {
        FALSE
      }
    })
    
    check_onclick_(NCleanData, 1, this_envir)
    
    
    output$error_maxN <- renderUI({
      if (is.na(input$maxN) || input$maxN > 0) {
        error_selection$error_maxN <- FALSE
        tags$span()
      }
      if (!is.na(input$maxN)) {
        if (input$maxN < 0) {
          error_selection$error_maxN <- TRUE
          tags$span("error: please select a value > 0",
            class = "errormsg")
        }
      }
    })
    
    
    lowComplexData <- reactive({
      if (!is.null(input$rm.lowcomplex)) {
        TRUE
      } else {
        FALSE
      }
    })
    
    check_onclick_(lowComplexData, 2, this_envir)
    
    output$error_complexThres <- renderUI({
      if (is.na(input$complexThres) || input$complexThres > 0) {
        error_selection$error_complexThres <- FALSE
        tags$span()
      }
      if (!is.na(input$complexThres)) {
        if (input$complexThres < 0) {
          error_selection$error_complexThres <- TRUE
          tags$span("error: please select a value > 0", 
            class = "errormsg")
        }
      }
    })
    
    
    SeqInput <- reactive({
      if (!is.null(input$rm.adapt)) {
        TRUE
      } else {
        FALSE
      }
    })
    
    output$error_errate <- renderUI({
      if (is.na(input$errate) || input$errate > 0) {
        error_selection$error_complexThres <- FALSE
        tags$span("")
      }
      if (!is.na(input$errate)) {
        if (input$errate < 0 || input$errate > 1) {
          error_selection$errate <- TRUE
          tags$span("error: please select a value between 0 and 1", 
            class = "errormsg")
        }
      }
    })
    
    check_onclick_(SeqInput, 3, this_envir)
    
    
    
    FilterbymeanQ <- reactive({
      if (!is.null(input$meanQ) && !is.null(input$thresmeanQ)) {
        TRUE
      } else {
        FALSE
      }
    })
    
    output$error_thresmeanQ <- renderUI({
      if (is.na(input$thresmeanQ) || input$thresmeanQ > 0) {
        error_selection$error_thresmeanQ <- FALSE
        tags$span("")
      }
      if (!is.na(input$thresmeanQ)) {
        if (input$thresmeanQ < 0) {
          error_selection$error_thresmeanQ <- TRUE
          tags$span("error: please select a value < 0", 
            class = "errormsg")
        }
      }
    })
    
    check_onclick_(FilterbymeanQ, 4, this_envir)
    
    
    
    TrimCleanData <- reactive({
      if (!is.null(input$rm.qual) && !is.null(input$thresQual)) {
        TRUE
      } else {
        FALSE
      }
    })
    
    
    output$error_thresQual <- renderUI({
      if (is.na(input$thresQual) || input$thresQual > 0) {
        error_selection$error_thresQual <- FALSE
        tags$span("")
      }
      if (!is.na(input$thresQual)) {
        if (input$thresQual < 0) {
          error_selection$error_thresQual <- TRUE
          tags$span("error: please select a value < 0", class = "errormsg")
        }
      }
    })
    
    check_onclick_(TrimCleanData, 5, this_envir)
    
    
    FixedCleanData <- reactive({
      if ((!is.null(input$rmFixed) && 
          !is.null(input$rm3)) || 
          (!is.null(input$rmFixed) && 
              !is.null(input$rm5))) {
        TRUE
      } else {
        FALSE
      }
    })
    
    
    
    output$error_rm35 <- renderUI({
      error_selection$error_rm35 <- ifelse(is.na(input$rm3), 
        FALSE, 
        ifelse(input$rm3 >  0, FALSE, TRUE)) || ifelse(is.na(input$rm5), 
          FALSE, ifelse(input$rm5 > 
              0, FALSE, TRUE))
      if (error_selection$error_rm35) {
        tags$span("error: please select values > 0", class = "errormsg")
      } else {
        tags$span("")
      }
    })
    
    
    check_onclick_(FixedCleanData, 6, this_envir)
    
    
    
    SizeCleanData <- reactive({
      if ((!is.null(input$rmSize) && !is.null(input$rmMin)) || 
          (!is.null(input$rmSize) && 
              !is.null(input$rmMax))) {
        TRUE
      } else {
        FALSE
      }
    })
    
    
    output$error_rmMinMax <- renderUI({
      error_selection$error_rmMinMax <- ifelse(is.na(input$rmMin), 
        FALSE, ifelse(input$rmMin > 0, FALSE, TRUE)) || 
        ifelse(is.na(input$rmMax), 
          FALSE, ifelse(input$rmMax >  0, FALSE, TRUE))
      if (error_selection$error_rmMinMax) {
        tags$span("error: please select values > 0", class = "errormsg")
      } else {
        tags$span("")
      }
    })
    
    
    check_onclick_(SizeCleanData, 7, this_envir)
    
    
    
    DuplicCleanData <- reactive({
      if (!is.null(input$rm.duplic)) {
        TRUE
      } else {
        FALSE
      }
    })
    
    check_onclick_(DuplicCleanData, 8, this_envir)
    
    
    output$selection <- renderUI({
      if (input$fileTypeIn == "SR") {
        actionButton(class = "setButton", "selection", "File")
      } else {
        tagList(actionButton(class = "setButton",
          "selection", "File 1"), 
          actionButton(class = "setButton", 
            "selection2", "File 2"))
      }
    })
    
    
    selection_detect <- reactiveValues(x = FALSE, y = FALSE)
    
    observeEvent(input$selection, {
      filepath$x <- try(file.choose(), silent = TRUE)
      selection_detect$x <- TRUE
      if (length(filepath$x) == 0) {
        filepath$x <- ""
        output$myPath <- renderPrint(cat("choose a file..."))
        output$encoding <- renderPrint(cat("..."))
        selection_detect$x <- FALSE
      }
    })
    
    observeEvent(input$selection2, {
      filepath$x[2] <- try(file.choose(), silent = TRUE)
      selection_detect$y <- TRUE
      if (!selection_detect$x) {
        filepath$x <- ""
        showModal(modalDialog(title = "FastqCleaner message",
          "Please load first File 1", 
          footer = modalButton("Ok")))
      }
    })
    
    my_encoding <- reactiveValues(value = NULL)
    
    observe({
      
      if (input$fileTypeIn == "SR") {
        req(selection_detect$x)
      } else {
        req(selection_detect$x)
        req(selection_detect$y)
      }
      
      isolate({
        progress <- Progress$new(session, min = 0, max = length(filepath$x))
        progress$set(message = "Analyzing file...")
        
        printPath <- TRUE
        printEncoding <- TRUE
        sampleQual <- integer()
        
        cat("Checking file and encoding...\n")
        t1 <- proc.time()
        
        if (length(filepath$x) > 1) {
          nsamp <- round(input$sample)/2
        } else {
          nsamp <- input$sample
        }
        
        for (i in filepath$x) {
          temporalPath <- try(ShortRead::FastqStreamer(i, n = nsamp), 
            silent = TRUE)
          thisTemporalFile <- try(ShortRead::yield(temporalPath), 
            silent = TRUE)
          try(close(temporalPath), silent = TRUE)
          
          if (length(thisTemporalFile) == 0 || 
              class(thisTemporalFile) == "try-error") {
            filepath$x <- ""
            showModal(modalDialog(title = "FastqCleaner message",
              "Invalid file: ShortRead cannot recognize the file format, 
              or your file is empty", 
              footer = modalButton("Ok")))
            printPath <- FALSE
            printEncoding <- FALSE
            break
          }
          
          new_sample <- Biostrings::quality(Biostrings::quality(thisTemporalFile))
          new_sample <- utf8ToInt(as.character(unlist(new_sample)))
          sampleQual <- c(sampleQual, new_sample)
          rm(thisTemporalFile)
        }
        my_encoding$value <- check_encoding(x = sampleQual)
        encodType <- "Encoding auto detected"
        
        t2 <- proc.time()
        outtime <- t2 - t1
        cat("OK...checks performed in ", outtime[3], " seconds\n\n")
        if (input$customEncod != 0) {
          my_encoding_custom <- check_encoding(custom = input$customEncod)
          if (!identical(my_encoding, my_encoding_custom)) {
            warning(paste0("The program detected the following encoding: ", 
              my_encoding, "\n", " You selected: my_encoding_custom"))
          }
          
          my_encoding$value <- my_encoding_custom
          encodType <- "Encoding selected by user"
        }
        
        radio_to_check <- my_encoding$value[["y"]]
        output$set_encoding_radio <- renderText({
          out <- list()
          for (i in 1:6) {
            if (radio_to_check != i) {
              out[[i]] <- HTML({
                paste0("<script>
                  $( \"input[name='customEncod'][value ='", 
                  i, "']\" ).prop(\"checked\", false);
                  </script>")
              })
              } else {
                out[[i]] <- HTML({
                  paste0("<script>
                    $( \"input[name='customEncod'][value ='", 
                    i, "']\" ).prop(\"checked\", true);
                    </script>")
                })
                }
            }
          unlist(out)
        })  
        progress$close()
        
        if (printPath) {
          output$myPath <- renderPrint(cat(length(filepath$x), 
            ifelse(length(filepath$x) == 1, 
              "file selected", 
              "files selected")))
        }
        
        
        if (printEncoding) {
          output$encoding <- renderPrint(cat("Encoding:", 
            my_encoding$value[["x"]], 
            "; sample range: ", my_encoding$value[["range"]], 
            "; n = ", input$nfile, 
            "\n", encodType))
        }
        
        selection_detect$x <- FALSE
        selection_detect$y <- FALSE
      })
      })
    
    observe({
      
      if (input$goButton == 0) {
        return()
      }
      isolate({
        
        if (any(filepath$x == "")) {
          
          showModal(modalDialog(title = "FastqCleaner message",
            "Please select a file first!", 
            footer = modalButton("Ok")))
          return()
        }
        
        errors_in_selection <- unlist(reactiveValuesToList(error_selection))
        if (any(errors_in_selection)) {
          showModal(modalDialog(title = "FastqCleaner message", 
            "Invalid values present. Please revise the entered values", 
            footer = modalButton("Ok")))
          return()
        }
        
        vars <- c("nFilterCatch", "lowComplexFilterCatch", "seqFilterCatch", 
          "meanQFilterCatch", "trim3FilterCatch",
          "fixedFilterCatch", "lengthFilterCatch", 
          "duplicFilterCatch")
        
        for (i in vars) {
          messages[[i]] <- c("Filter unused", "")
          readsWidth[[i]] <- list(list(numeric(), numeric()), list(numeric(), 
            numeric()))
        }
        messages$outResult = c("No file processed yet", "")
        
        
        proctime_start <- proc.time()
        
        
        processingFunction_(this_envir)
        
        
        proctime_end <- proc.time()
        
        proctime <- proctime_end - proctime_start
        
        cat("-------------------------------\n")
        cat("\nProcessing finished\n\n")
        cat("Processing performed in ", proctime[1], "seconds\n")
        
        
        cat("Processing successful!\n")
        showModal(modalDialog(title = "FastqCleaner message", "Done!", 
          footer = modalButton("Ok")))
        hasRun$x <- TRUE
      })
      
    })  
    
    reactP2 <- reactive(reactiveValues(filepath = filepath$x,
      meanQ = input$meanQ,
      thresmeanQ = input$thresmeanQ, 
      rm.qual = input$rm.qual, 
      thresQual = input$thresQual,
      rmFixed = input$rmFixed, 
      rm3 = input$rm3, 
      rm5 = input$rm5, 
      rMin = input$rmMin, 
      rmMax = input$rmMax))
    
    
    callModule(panel1, "panel1", messages = reactive(messages), rList = reactP2)
    
    
    reactP3 <- reactive(reactiveValues(filepath = filepath$x,
      fileTypeIn = input$fileTypeIn, 
      goButton = input$goButton))
    
    
    temporalPlot <- reactiveValues(x = list())
    
    isPaired <- reactiveValues()
    observe({
      if (input$fileTypeIn == "SR") {
        isPaired$x <- FALSE
      } else {
        isPaired$x <- TRUE
      }
    })
    
    
    observeEvent(input$fileTypeIn, 
      temporalPlot$x <- list(HTML("<script></script>")))
    
    error_state <- reactiveValues(x = FALSE)
    
    error_state_init <- reactiveValues(x = FALSE)
    
    empty_plot <- reactiveValues(x = FALSE)
    
    has_new_plots <- reactiveValues(x = FALSE)
    
    validate_plot <- debounce(reactive({
      hasRun$x
      input$filePlot
      input$sample
      input$klength
      input$maxFreq
      
      has_new_plots$x <- FALSE
      empty_plot$x <- FALSE
      
      updateSelectInput(session, "plotType", selected = "...")
      
      
      if (input$goButton == 0 || any(filepath$x == "")) {
        temporalPlot$x <- list(HTML("<script></script>"))
        return(FALSE)
      }
      
      error_state$x <- !(FastqCleaner:::isNaturalNumber(list(input$sample - 1, 
        input$klength, input$maxFreq)))
      
      if (error_state$x) {
        updateSelectInput(session, "plotType", selected = "...")
        output$plot1 <- renderPrint(HTML("<script>
          $( '#container1' ).addClass('alert');
          $( '#container1' ).append( '<p> Oops... sample size, k-mer length and top sequences must be integer numbers! </p>' );
          $( '.cssload-fond' ).addClass('hidden');
          </script>"))
        error_state_init$x <- TRUE  
        return(FALSE)
        
      }
      
      if (error_state_init$x) {
        output$plot1 <- renderPrint(HTML("<script>
          $( '#container1' ).removeClass('alert');
          $( '#container1' ).children().remove();
          $( '.cssload-fond' ).removeClass('hidden');
          </script>"))
        error_state_init$x <- FALSE
      }
      
      
      if (input$filePlot == "my_output") {
        if (!isPaired$x) {
          length_try <- length(ShortRead::yield(
            temp_stream <- try(
              ShortRead::FastqStreamer(paste0(filepath$x[1], 
                "_out"), 10)), silent = TRUE))
          if (length_try == 0) {
            try(close(temp_stream), silent = TRUE)
            empty_plot$x <- TRUE
            return(TRUE)
          }
        } else {
          length_try_L <- length(ShortRead::yield(
            temp_stream <- try(
              ShortRead::FastqStreamer(paste0(filepath$x[1], 
                "_out"), 10)), silent = TRUE))
          if (length_try_L == 0) {
            try(close(temp_stream), silent = TRUE)
            empty_plot$x <- TRUE
            return(TRUE)
          }
          length_try_R <- length(ShortRead::yield(
            temp_stream <- try(
              ShortRead::FastqStreamer(paste0(filepath$x[1], 
                "_out"), 10)), silent = TRUE))
          if (length_try_R == 0) {
            try(close(temp_stream), silent = TRUE)
            empty_plot$x <- TRUE
            return(TRUE)
          }
        }
      }
      
      
      output$plot1 <- renderPrint(HTML("<script>
        $( '#container1' ).children().remove();
        </script>"))
      
      TRUE
    }), 1000)
    
    
    observe({
      validate_plot
      if (!validate_plot()) {
        return()
      }
      
      
      isolate({
        updateSelectInput(session, "plotType", selected = "...")
        
        progress <- Progress$new(session)
        progress$set(message = "Computing plots...")
        
        if (input$filePlot == "my_input") {
          temporalPlot$x <- FastqCleaner:::myPlot(isPaired$x, 
            filepath$x, input$sample, 
            input$klength, "input", input$maxFreq)
          
        } else if (input$filePlot == "my_output") {
          
          temporalPlot$x <- FastqCleaner:::myPlot(isPaired$x, 
            sprintf("%s_out",filepath$x), 
            input$sample, 
            input$klength,
            "output", input$maxFreq)
          
        }
        progress$close()
        has_new_plots$x <- TRUE
      })
      
    })
    
    
    output$plotPanel <- renderUI({
      
      
      cuantos <- input$fileTypeIn
      if (cuantos == "SR") {
        nplots <- 1
      } else {
        nplots <- 2
      }
      
      
      if (!hasRun$x || input$goButton == 0 || 
          any(filepath$x == "") || error_state$x || 
          empty_plot$x) {
        outplots <- tags$div(id = "container1", 
          style = "height: 400px; margin: auto; width: 800px")
        output$tabBut <- renderUI(HTML("<p></p>"))
      } else {
        
        outplots <- lapply(1:nplots, function(i) {
          plotname <- paste0("plot", i)
          tags$div(id = paste0("container", i), 
            style = "height: 400px; margin: auto; width: 800px")
        })
        outplots <- do.call(tagList, outplots)
      }
      
      div(id = "plot-container", 
        style = "margin-left:150px", 
        div(class = "cssload-fond", 
          div(class = "cssload-container-general", 
            div(class = "cssload-internal", 
              div(class = "cssload-ballcolor cssload-ball_1")), 
            div(class = "cssload-internal", 
              div(class = "cssload-ballcolor cssload-ball_2")), 
            div(class = "cssload-internal", 
              div(class = "cssload-ballcolor cssload-ball_3")),
            div(class = "cssload-internal", 
              div(class = "cssload-ballcolor cssload-ball_4")))), 
        outplots)
      
    })
    
    
    observeEvent(input$plotType, {
      
      if (!has_new_plots$x) {
        return()
      }
      
      if (input$goButton == 0 || any(filepath$x == "") || error_state$x) {
        
        if (!error_state$x) {
          output$plot1 <- renderPrint(HTML("<script></script>"))
        }
        output$plot2 <- renderPrint(HTML("<script></script>"))
        output$tabBut <- renderUI(HTML("<script></script>"))
        output$freqTable1 <- DT::renderDataTable(data.frame())
        output$freqTable2 <- DT::renderDataTable(data.frame())
        return()
      }
      
      isolate({
        
        if (empty_plot$x) {
          
          output$plot1 <- renderPrint(HTML("<script>
            $( '#container1' ).addClass('alert');
            $( '.highcharts-container ' ).empty();
            $( '#container1' ).append( '<p>Oops... the output file is empty! </p>' );
            $( '.cssload-fond' ).addClass('hidden');
            </script>"))
          updateSelectInput(session, "plotType", selected = "...")
          
          return()
        }
        
        if (is.null(input$plotType)) {
          selected <- "a"
        } else {
          selected <- input$plotType
        }
        
        
        output$plot1 <- renderPrint(tagList(HTML("<script>
          $('#container1').removeClass('alert');
          $('#container1').children().remove();
          $('.cssload-fond').removeClass('hidden');
          </script>"), 
          HTML(temporalPlot[["x"]][[1]][[selected]])))
        
        if (input$fileTypeIn == "SR") {
          
          
          if (selected == "j") {
            output$freqTable1 <- DT::renderDataTable(
              temporalPlot$x[[1]][["oligofreq"]], 
              width = "20%", 
              options = list(lengthChange = FALSE, pageLength = 8))
            output$tabBut <- renderUI(
              actionButton("tabBut", "View kmer table"))
          } else if (selected == "i") {
            output$freqTable1 <- DT::renderDataTable(
              temporalPlot$x[[1]][["table"]], 
              width = "20%", options = list(lengthChange = FALSE, 
                pageLength = 8))
            output$tabBut <- renderUI(actionButton("tabBut", "View Table"))
          } else {
            
            output$tabBut <- renderUI(HTML("<p></p>"))
          }
          
        } else {
          
          output$plot1 <- renderPrint(HTML(temporalPlot[["x"]][[1]][[selected]]))
          output$plot2 <- renderPrint(HTML(temporalPlot[["x"]][[2]][[selected]]))
          if (selected == "j") {
            output$freqTable1 <- DT::renderDataTable(
              temporalPlot[["x"]][[1]][["oligofreq"]], 
              width = "20%", options = list(lengthChange = FALSE, 
                pageLength = 8))
            output$freqTable2 <- DT::renderDataTable(
              temporalPlot[["x"]][[2]][["oligofreq"]], 
              width = "20%", options = list(lengthChange = FALSE, 
                pageLength = 8))
            output$tabBut <- renderUI(actionButton("tabBut", 
              "View kmer table"))
          } else if (selected == "i") {
            output$freqTable1 <- DT::renderDataTable(
              temporalPlot[["x"]][[1]][["table"]], 
              width = "20%", options = list(lengthChange = FALSE, 
                pageLength = 8))
            output$freqTable2 <- DT::renderDataTable(
              temporalPlot[["x"]][[2]][["table"]], 
              width = "20%", options = list(lengthChange = FALSE, 
                pageLength = 8))
            output$tabBut <- renderUI(actionButton("tabBut", "View Table"))
          } else {
            output$tabBut <- renderUI(HTML("<p></p>"))
          }
          
        }
      })
    })
    
    callModule(panel2, "panel2")
    
    callModule(panel3, "panel3", arg = reactive(input$navbar))
    
    })

