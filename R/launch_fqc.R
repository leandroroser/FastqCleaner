#' Launch FastqCleaner application
#' @param launch.browser Launch in browser? Default TRUE
#' @param ... Additional parameters passed to \code{\link[shiny]{runApp}}
#' @author Leandro Roser \email{learoser@@gmail.com}
#' @examples 
#' # Paste in te console to launch the application:
#' # launch_fqc() 
#' 
#' NULL
#' @return Launch the application, without return value
#' @export


launch_fqc <- function(launch.browser = TRUE, ...) {
    fpath <- system.file("application", package = "FastqCleaner")
    shiny::runApp(fpath, launch.browser = launch.browser, ...)
}
