#' Plot mean absolute SHAP values for a gene
#' @param results A `SignatureResults` object.
#' @param gene Gene to plot.
#' @param comparison_id Optional comparison ID.
#' @param top_n Number of features to show.
#' @return A ggplot object.
#' @export
plot_shap <- function(results,gene,comparison_id=NULL,top_n=20L) {
    if (!inherits(results,"SignatureResults")) stop("results must be a SignatureResults object.")
    dat <- results$shap_summary[results$shap_summary$gene==gene,,drop=FALSE]
    if (!nrow(dat)) stop("No SHAP results found for gene: ",gene)
    if (!is.null(comparison_id)) dat <- dat[dat$comparison_id==comparison_id,,drop=FALSE]
    else if (length(unique(dat$comparison_id))>1L) stop("Multiple comparisons found for ",gene,". Specify comparison_id.")
    if (!nrow(dat)) stop("No SHAP results found for the requested comparison.")
    dat <- head(dat[order(dat$mean_abs_shap,decreasing=TRUE),,drop=FALSE],top_n); dat$feature <- factor(dat$feature,levels=rev(dat$feature))
    ggplot2::ggplot(dat,ggplot2::aes(x=feature,y=mean_abs_shap)) + ggplot2::geom_col() + ggplot2::coord_flip() +
        ggplot2::labs(title=paste(gene,unique(dat$positive),"vs",unique(dat$negative)),x=NULL,y="Mean |SHAP|") + ggplot2::theme_bw()
}
