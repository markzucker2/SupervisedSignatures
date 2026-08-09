#' Create an empty SignatureResults object
#' @param metadata Optional named list of analysis metadata.
#' @return An object of class `SignatureResults`.
#' @export
SignatureResults <- function(metadata = list()) {
    out <- list(summary=data.frame(), models=list(), predictions=list(),
                shap_raw=list(), shap_summary=data.frame(), metadata=metadata)
    class(out) <- "SignatureResults"
    out
}

#' @export
print.SignatureResults <- function(x, ...) {
    cat("<SignatureResults>\n")
    cat("Models:", length(x$models), "\n")
    cat("Summary rows:", nrow(x$summary), "\n")
    cat("SHAP summary rows:", nrow(x$shap_summary), "\n")
    invisible(x)
}
