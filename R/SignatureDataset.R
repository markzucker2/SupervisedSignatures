#' Create a SignatureDataset object
#'
#' @param samples Data frame containing sample identifiers. Must contain `sample_id`.
#' @param features Data frame of model features. Must contain `sample_id`;
#'   every other column is treated as a feature.
#' @param covariates Optional data frame of covariates. Must contain `sample_id`.
#' @param gene_states Data frame containing gene states. Must contain
#'   `sample_id`, `gene`, and `state`.
#' @param feature_type Character label describing the feature representation.
#' @param feature_set Character label describing the feature set.
#' @param metadata Optional named list of additional metadata.
#' @return An object of class `SignatureDataset`.
#' @export
SignatureDataset <- function(samples, features, covariates = NULL,
                              gene_states, feature_type = "unspecified",
                              feature_set = "unspecified", metadata = list()) {
    validate_required_columns(samples, "sample_id", "samples")
    validate_required_columns(features, "sample_id", "features")
    validate_required_columns(gene_states, c("sample_id", "gene", "state"), "gene_states")
    if (!is.null(covariates)) validate_required_columns(covariates, "sample_id", "covariates")

    sample_ids <- unique(as.character(samples$sample_id))
    if (!all(unique(as.character(features$sample_id)) %in% sample_ids))
        stop("All feature sample_id values must occur in samples.")
    if (!all(unique(as.character(gene_states$sample_id)) %in% sample_ids))
        stop("All gene_states sample_id values must occur in samples.")
    if (!is.null(covariates) && !all(unique(as.character(covariates$sample_id)) %in% sample_ids))
        stop("All covariate sample_id values must occur in samples.")

    out <- list(
        samples = samples,
        features = features,
        covariates = covariates,
        gene_states = gene_states,
        metadata = c(list(feature_type = feature_type, feature_set = feature_set), metadata)
    )
    class(out) <- "SignatureDataset"
    out
}

#' @export
print.SignatureDataset <- function(x, ...) {
    cat("<SignatureDataset>\n")
    cat("Samples:", nrow(x$samples), "\n")
    cat("Features:", ncol(x$features) - 1L, "\n")
    cat("Genes:", length(unique(x$gene_states$gene)), "\n")
    cat("Feature type:", x$metadata$feature_type, "\n")
    cat("Feature set:", x$metadata$feature_set, "\n")
    invisible(x)
}
