
  
 navbarPage(
    id = "navbar",
    title = list(tags$div(
      tags$img(src = "images/icon.png",
        alt = "icon"),
      tags$p(" FastqCleaner")
    )),
    header = tags$head(
          tags$link(
              rel = "icon",
              href = "images/icon.png",
              type = "image/png",
              width = "1em",
              height = "1em"
          ),
          tags$title("FastqCleaner"),
          tags$link(rel = "stylesheet",
              type = "text/css",
              href = "styles/masterCss.css")
      ),
    
    
    tabPanel( "Clean the Data", sidebarLayout(
      
      sidebarPanel(
        width = 3,
        class = "sideBar",
        h3("FastqCleaner"),
        p("A program to clean FASTQ files"),
        br(),
        hr(),
        br(),
        div() 
        
      ),
      
      mainPanel(
        
        fluidRow(
          br(),
          
          column(6, class = "left",
            singleton(tags$head(
              tags$script(
                'Shiny.addCustomMessageHandler("testmessage",
                function(message) { alert(message);
                });'                               )
      )),
            h3(class = "blackTag", "Select the operations")
            ),
          
          column(3, class = "right"),
          br(),
          br(),
          br(),
          br()
            ), 
        
        bsPopover("help", "",
          "Help - documentation",
          options = list(container = "body")
        ),
        
        
        fluidRow(
          
          column(6, class = "accordion-container",
            tags$button(id = "button1", class = "accordion", 
              "1. Remove by N(s) #",
              tags$span()
            ),
            
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rm.N",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                numericInput("maxN",
                  label = p(class = "controlcheck-p1", 
                    "Maximum number of N(s)"),
                  value = 3,
                  min = 0
                ),
                uiOutput("error_maxN")
              )
            ), 
            
            bsPopover("maxN", "",
              "The maximum number of Ns in a read",
              options = list(container = "body")
            ),
            
            
            tags$button(
              id = "button2",
              class = "accordion",
              "2. Remove low complexity sequences",
              tags$span()
            ),
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rm.lowcomplex",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                numericInput("complexThres",
                  label = p(class = "controlcheck-p1", 
                    "Complexity threshold"),
                  value = 0.5,
                  min = 0
                ),
                uiOutput("error_complexThres")
              )
            ),
            
            bsPopover("complexThres","",
              "Shanon entropy/reference entropy",
              options = list(container = "body")
            ),
            
            
            tags$button(id = "button3", class = "accordion",
              "3. Remove adapters",
              tags$span()
            ),
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rm.adapt",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                
                textInput("LpatternF",
                  label = p(class = "controlcheck-p1",
                    "Left forward adapter"),
                  value = ""
                ),
                checkboxInput("reverseLF",
                  label = p(class = "controlcheck-p1-radio",
                    "use reverse complement?"),
                  value = FALSE
                ),
                textInput("RpatternF",
                  label = p(class = "controlcheck-p1",
                    "Right forward adapter"),
                  value = ""
                ),
                checkboxInput("reverseRF",
                  label = p(class = "controlcheck-p1-radio",
                    "use reverse complement?"),
                  value = FALSE
                ),
                textInput("LpatternR",
                  label = p(class = "controlcheck-p1",
                    "Left reverse adapter"),
                  value = ""
                ),
                checkboxInput("reverseLR",
                  label = p(class = "controlcheck-p1-radio",
                    "use reverse complement?"),
                  value = FALSE
                ),
                textInput("RpatternR",
                  label = p(class = "controlcheck-p1",
                    "Right reverse adapter"),
                  value = ""
                ),
                checkboxInput("reverseRR",
                  label = p(class = "controlcheck-p1-radio",
                    "use reverse complement?"),
                  value = FALSE
                ),
                
                
                tags$button(id = "advancedAdapterButton", 
                  class = "adButton",
                  "Advanced",
                  style="margin-left:2em"
                ),
                br(), br(),
                
                
                fluidRow(id = "advancedAdapterMenu",
                  column(12,style = "float:left;",
                    checkboxInput("indels",
                      label = p(class = "controlcheck-p1-radio",
                        "with Indels?"),
                      value = FALSE
                    ),
                    checkboxInput("anchored",
                      label = p(class = "controlcheck-p1-radio",
                        "anchored adapters?"),
                      value = TRUE
                    ),
                    
                    numericInput("min_match_flank",
                      label = p(class = "controlcheck-p1",
                        "Minimum match in flanks"),
                      value = 1,
                      min = 0,
                      step = 1
                    ),
                    
                    
                    numericInput("errate",
                      label = p(class = "controlcheck-p1",
                        "Error rate"),
                      value = 0.2,
                      min = 0,
                      max = 1,
                      step = 0.1
                    ),
                    uiOutput("error_errate"),
                    
                    
                    radioButtons("first_forward",
                      label = p(class = "controlcheck-p1", 
                        "Remove first in PE forward"),
                      choices = list("5'" = "L", "3'" = "R"),
                      selected = "R"
                    ),
                    radioButtons("first_reverse",
                      label = p(class = "controlcheck-p1", 
                        "Remove first in PE reverse"),
                      choices = list("5'" = "L", "3'" = "R"),
                      selected = "R"
                    )
                  ) 
                )
              )
              
            ),
            
            bsPopover("LpatternF", "",
              "The left pattern for forward sequences",
              options = list(container = "body")
            ),
            
            
            bsPopover("RpatternF", "",
              "The right pattern for forward sequences",
              options = list(container = "body")
            ),
            
            bsPopover("LpatternR", "",
              "The left pattern for reverse sequences",
              options = list(container = "body")
            ),
            
            bsPopover("RpatternR", "",
              "The right pattern for reverse sequences",
              options = list(container = "body")
            ),
            
            bsPopover("indels", "", 
              "allow indels?",
              options = list(container = "body")
            ),
            
            bsPopover("anchored", "", 
              "anchored adapters?",
              options = list(container = "body")
            ),
            
            bsPopover("min_match_flank", "",
              "minimum overlap between adapter and reads 
              to trim for a match in flanking region of reads",
              options = list(container = "body")
            ),
            
            bsPopover("errate", "",
              "error rate value for error rate method",
              options = list(container = "body")
            ),
            
            
            tags$button(id = "button4", class = "accordion",
              "4. Filter by average quality?",
              tags$span()
            ),
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("meanQ",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                
                numericInput("thresmeanQ",
                  label = p(class = "controlcheck-p1",
                    "Average Quality threshold"),
                  value = 20,
                  min = 0
                ),
                uiOutput("error_thresmeanQ")
              )
              
            ),
            
            bsPopover("thresmeanQ", "",
              "Mean quality threshold",
              options = list(container = "body")
            )
            
          ), 
          
          column(6, class = "accordion-container",
            
            
            tags$button(id = "button5", class = "accordion",
              "5. Trim low quality 3' tails?",
              tags$span()
            ),
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rm.qual", label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                numericInput("thresQual",
                  label = p(class = "controlcheck-p1", 
                    "Quality threshold"),
                  value = 20,
                  min = 0
                ),
                uiOutput("error_thresQual")
                
              )
            ), 
            
            bsPopover("thresQual", "",
              "Threshold quality in  tail",
              options = list(container = "body")
            ),
            
            tags$button(id = "button6", class = "accordion",
              "6. Trim 3' or 5' by a fixed number?",
              tags$span()
            ),
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rmFixed",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                numericInput("rm3",
                  label = p(class = "controlcheck-p1", 
                    "Trim 3'"),
                  value = NULL,
                  min = 0
                ),
                numericInput("rm5",
                  label = p(class = "controlcheck-p1", 
                    "Trim 5'"),
                  value = NULL,
                  min = 0
                ),
                uiOutput("error_rm35")
              )
            ), 
            
            bsPopover("rm3", "",
              "Number of bases to trim from 3 end",
              options = list(container = "body")
            ),
            
            bsPopover("rm5", "",
              "Number of bases to trim from 5 end",
              options = list(container = "body")
            ),
            
            tags$button(id = "button7", class = "accordion",
              "7. Filter sequences by length?",
              tags$span()
            ),
            
            div(class = "panelCustom",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rmSize",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                ),
                numericInput("rmMin",
                  label = p(class = "controlcheck-p1", 
                    "Min length"),
                  value = NULL,
                  min = 0
                ),
                numericInput("rmMax",
                  label = p(class = "controlcheck-p1", 
                    "Max length"),
                  value = NULL,
                  min = 0
                ),
                uiOutput("error_rmMinMax")
              )
            ), 
            
            bsPopover("rmMin", "",
              "Minimum number of bases in a read",
              options = list(container = "body")
            ),
            
            bsPopover("rmMax", "",
              "Maximum number of bases in a read",
              options = list(container = "body")
            ),
            
            tags$button(id = "button8", class = "accordion",
              "8. Remove duplicated sequences?",
              tags$span()
            ),
            
            div(class = "panelCustom",
              style = "margin-bottom:20em",
              wellPanel(class = "inAccordion",
                checkboxGroupInput("rm.duplic",
                  label = "",
                  choices = list("Use filter?" = 1),
                  selected = 0
                )
              )
            ) 
          ),
          
          
          br(), br()
          
          ), 
        
        
        br(), br(), br(),
        
        
        fluidRow(
          column(6, class = "left",
            h3(class = "blackTag", "Select a file")),
          column(3, class = "right"),
          br(),br(),br(),br()
        ),
        
        
        fluidRow(id = "filecontainer",
          column(6,
            fluidRow(
              column(12, style = "float:left; margin-left:2em;margin-top:1em;",
                radioButtons("fileTypeIn", 
                  p(class="controlcheck-p3", "Library type"),
                  choices = c("single-end reads" = "SR", 
                    "paired-end reads" = "PE"))
              ),
              column(12, style = "float:left; margin-left:2em;",
                htmlOutput("selection"),
                actionButton(class = "setButton", "goButton", "Run!")
              ),
              
              column(12, style = "float:left; margin-left:2em; margin-top: 2em",
                p(class="controlcheck-p3", "File selected:"),
                verbatimTextOutput("myPath")
              ),
              column(12, style = "float:left; margin-left:2em;margin-top:1em;",
                p(class="controlcheck-p3", "Encoding:"),
                verbatimTextOutput("encoding")
              )
            )
          ), 
          
          
          column(6,
            fluidRow(
              column(12, style = "float:left; margin-left:2em;",
                style = "float:left; margin-right:2em;;margin-top:1em;",
                radioButtons("fileTypeOut",
                  p(class="controlcheck-p3", "Output format"),
                  choices = c("FASTQ" = "fastq", "gz (compressed)" = "gz"),
                  selected = "fastq")
              ),
              column(12, style = "float:left; margin-left:2em;",
                actionButton(class = "setButton", "resetSelections", "Clear"),
                actionButton(class = "setButton", "resetAll", "Reset")
              )
            )
          ),
          
          
          
          fluidRow(
            column(12, 
              tags$button(id = "advancedButton",
                class = "adButton",
                "Advanced"
              )
            ),
            
            column(12,
              fluidRow(id = "advancedMenu",
                column(6, style = "float:left; margin-left:2em;",
                  radioButtons("customEncod",
                    p(class ="controlcheck-p1 encoding",
                      "The program will auto detect encoding.
                      You can also select a custom format of the following:"),
                    choices = c(
                      "None selected" = 0,
                      "Sanger    | Q range: [33-73]" = 1,
                      "Illumina 1.8+   | Q range: [33-74]" = 2,
                      "Illumina 1.3+   | Q range: [64-104]" = 3,
                      "Illumina 1.5+   | Q range: [66-104]" = 4,
                      "Solexa    | range: [59-104]" = 5
                    )
                    ) 
              ), 
                column(6, style = "float:left; margin-left:2em;",
                  numericInput("nfile",
                    label = p(class ="controlcheck-p1 nfile", 
                      "Reads sampled for encoding detection and file processing"),
                    value = 100000
                  )
                )
                
              ) 
            )
          ) 
        ), 
        
        
        
        tags$script(src = 'scripts/javascript.js'),
        
        singleton(tags$head(
          tags$script(src = 'shared/jqueryui/jquery-ui.min.js')
        )),
        
        tags$script(src = "scripts/highcharts.js"),
        tags$script(src = "scripts/highcharts-more.js"),
        tags$script(src = "scripts/exporting.js"),
        tags$script(src = "scripts/offline-exporting.js"),
        tags$script(src = "scripts/FqC_highcharts_style.js"),
        tags$script(src = "scripts/heatmap.js"),
        
        htmlOutput('archivo'),
        htmlOutput('check1'),
        htmlOutput('check2'),
        htmlOutput('check3'),
        htmlOutput('check4'),
        htmlOutput('check5'),
        htmlOutput('check6'),
        htmlOutput('check7'),
        htmlOutput('check8'),
        htmlOutput('runDone'),
        htmlOutput('set_encoding_radio')
        
            ) 
    ) 
      ), 
    
    panel1Input("panel1"),
    
    
    tabPanel("Live Results", sidebarLayout(
      
      sidebarPanel(width = 3, id = "menuPlot", class = "sideBarPlot",
        numericInput("sample", label = p(class = "controlcheck-p3", 
          "Sample size"),  value = 10000),
        radioButtons("filePlot",
          p(class = "controlcheck-p3", "File to plot:"),
          choices = c("Input" = "my_input", "Output" = "my_output")
        ),
        br(),
        selectInput("plotType",
          label = p(class = "controlcheck-p3", "Diagnostic plots"),
          choices = c(
            "Select a plot..." = "",
            "Per cycle quality" = "a",
            "Per cycle mean base quality" = "b",
            "Mean quality distribution" = "c",
            "% reads with Phred scores > threshold" = "d",
            "Per cycle base proportion (barplot)" = "e",
            "Per cycle base proportion (lineplot)" = "f",
            "CG content over all reads" = "g",
            "Read length" = "h",
            "Read occurrence" = "i",
            "Relative k-mer diversity" = "j"
            
            
          ),
          selected = "none"
        ), 
        br(),
        numericInput(
          "klength",
          label = p(class = "controlcheck-p3", "Select kmer size"),
          value = 6
        ),
        br(),
        numericInput(
          "maxFreq",
          label = p(class = "controlcheck-p3", 
            "Top sequences in duplication level analysis"),
          value = 20
        ),
        
        br(),br(),br(),br(),br(),br(),br(),br(),
        br(),br(),br(),br(),br(),br(),br(),br()
      ),
      
      
      mainPanel(
        htmlOutput("plot1"),
        htmlOutput("plot2"),
        uiOutput("plotPanel"),
        bsModal("modalTable",
          title = "Most represented sequences",
          trigger = "tabBut",
          size = "large",
          DT::dataTableOutput('freqTable1'),
          DT::dataTableOutput('freqTable2')
        ),
        uiOutput("tabBut")
      )
      
    ) 
    ), 
    
    panel2Input("panel2"),
    
    panel3Input("panel3")
    
) 