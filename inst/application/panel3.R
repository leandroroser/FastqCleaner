panel3 <- function(input, output, session, arg) {
    observe({
        if (arg() == "Stop") {
            stopApp()
        }
    })
    
}
