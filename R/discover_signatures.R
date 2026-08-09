#' Run supervised signature discovery across comparisons
#' @param dataset A `SignatureDataset`.
#' @param comparisons Optional comparison data frame. If `NULL`, defaults are generated with `make_gene_state_comparisons()`.
#' @param reference_state Reference state for automatic comparison generation.
#' @param states Optional altered states to include.
#' @param compute_shap Logical; compute SHAP for each fit.
#' @param ... Arguments passed to `fit_classifier()`.
#' @return A `SignatureResults` object.
#' @export
discover_signatures <- function(dataset,comparisons=NULL,reference_state="WT",states=NULL,compute_shap=TRUE,...) {
    if (!inherits(dataset,"SignatureDataset")) stop("dataset must be a SignatureDataset.")
    if (is.null(comparisons)) comparisons <- make_gene_state_comparisons(dataset$gene_states,reference_state,states)
    if (!all(c("gene","positive","negative") %in% names(comparisons))) stop("comparisons must contain gene, positive, and negative.")
    results <- SignatureResults(metadata=c(dataset$metadata,list(reference_state=reference_state,n_comparisons=nrow(comparisons))))
    if (!nrow(comparisons)) { warning("No comparisons were generated."); return(results) }
    for (i in seq_len(nrow(comparisons))) {
        cmp <- comparisons[i,,drop=FALSE]; analysis <- build_dataset(dataset,cmp); fit <- fit_classifier(analysis,...)
        if (isTRUE(compute_shap)) fit <- compute_shap(fit)
        key <- if ("id" %in% names(cmp)) as.character(cmp$id[[1L]]) else paste(analysis$gene,analysis$positive,"vs",analysis$negative,sep="_")
        results$models[[key]] <- fit$model; results$predictions[[key]] <- fit$predictions; results$summary <- rbind(results$summary,fit$metrics)
        if (!is.null(fit$shap_summary)) results$shap_summary <- rbind(results$shap_summary,fit$shap_summary)
        if (is.null(results$metadata$fits)) results$metadata$fits <- list(); results$metadata$fits[[key]] <- fit
    }
    results$metadata$comparisons <- comparisons; results
}
