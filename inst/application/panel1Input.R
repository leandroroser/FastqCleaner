
panel1Input <-  function (id) {
 ns <- NS(id)
 
 
 tabPanel("File Operations", sidebarLayout(
  
  sidebarPanel(width = 3, class = "sideBar",
               h3("FastqCleaner"),
               p("A program to clean FASTQ files"),
               br(),
               hr(),
               br(),
               div() # .well.sideBar div
  ),
  
  mainPanel(
   
   # input/output ---------------------------------
   fluidRow(
    br(),
    column(6, class = "left",
           h3(class = "blackTag", "Input/Output")),
    column(3, class = "right"),
    br(), br(), br(), br(), br()
   ),
   
   fluidRow(
    column(6,
           p(class = "results", "File input:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("filepath")))
    ),
    column(6,
           p(class = "results", "File output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("outResult")))
    ),
    br(), br(), br(), br(), br()
   ),
   
   # file operations ------------------------------
   fluidRow(
    br(),
    # esta funcion permite enviar mensajes
    column(6, class = "left",
           h3(class = "blackTag", 
              "File operations")
    ),
    column(3, class = "right"),
    br(), br(), br(), br(), br()
   ),
   
   # Outputs right ---------------------------------------------
   fluidRow(
    
    column(6,
           p(class = "results", "N filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("N"))),
           p(class = "results", "Low complexity filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("lowComplexFilter"))),
           p(class = "results", "Adapter filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("adaptFilter"))),
           p(class = "results", "Average quality filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("meanQFilter")))
    ),
    
    # outputs left --------------------------------------------------
    column(6,
           p(class = "results", "Low quality 3' tails filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("TrimQual"))),
           p(class = "results", "Fixed-length tails trim output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("trimTails"))),
           p(class = "results", "Length filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("sizeFilter"))),
           p(class = "results", "Duplicated filter output:"),
           div(style = "width:80%",
               verbatimTextOutput(ns("DuplicFilter")))
           )
    
   ) # end fluidRow
   
  ) # end mainpanel
 ) # end sidebarLayout
 ) # end tabpanel
}
