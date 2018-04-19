
panel2Input <-  function(id) {
 ns <- NS(id)

tabPanel("About", sidebarLayout(
 sidebarPanel(
  width = 3,
  class = "sideBar",
  h3("FastqCleaner"),
  p("A program to clean FASTQ files"),
  br(),
  hr(),
  div()
 ),
 
 mainPanel(
  br(),
  br(),
  tags$a(
   href = "https://github.com/leandroroser",
   target = "_blank",
   fluidRow(
    id = ns("leandroUp"),
    column(3, id = ns("about"),
           div(class = "textname", "Follow the project in GitHub!")),
    tags$img(id = ns("github"), src = "images/github.svg")
   )
  ),
  
  fluidRow(class = "leandroDown",
           br(),
           column(
            3,
            id = ns("myMail"),
            div("Feedback to:", style = "color:blue"),
            tags$a("learoser@gmail.com", href =
                    "mailto:learoser@gmail.com",
                   style = "color:black;font-weight:bold;")
           ))
 ) # end mainPanel
) # end sidebarLayout
) # end tabPanel
}
