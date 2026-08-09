#' Cluster genes by SHAP-derived signatures
#' @param results A `SignatureResults` object.
#' @param genes Optional genes to include.
#' @param comparison Optional comparison ID.
#' @param method Hierarchical clustering method.
#' @param distance Distance method.
#' @param k Optional number of clusters.
#' @return A `SignatureClusters` object.
#' @export
cluster_signatures <- function(results,genes=NULL,comparison=NULL,method="complete",distance="euclidean",k=NULL) {
    if (!inherits(results,"SignatureResults")) stop("results must be a SignatureResults object.")
    dat <- results$shap_summary; if (!is.null(genes)) dat <- dat[dat$gene %in% genes,,drop=FALSE]; if (!is.null(comparison)) dat <- dat[dat$comparison_id==comparison,,drop=FALSE]
    if (!nrow(dat)) stop("No SHAP summary data available.")
    gene_counts <- tapply(dat$comparison_id,dat$gene,function(x) length(unique(x))); if (any(gene_counts>1L)) stop("Multiple comparisons exist for at least one gene. Specify comparison.")
    mat <- reshape_to_signature_matrix(dat); d <- stats::dist(mat,method=distance); hc <- stats::hclust(d,method=method); clusters <- NULL
    if (!is.null(k)) { if (k<2L || k>nrow(mat)) stop("k must be between 2 and the number of genes."); clusters <- stats::cutree(hc,k=k) }
    out <- list(signature_matrix=mat,distance=d,hclust=hc,clusters=clusters); class(out) <- "SignatureClusters"; out
}

#' @export
plot.SignatureClusters <- function(x,...) { graphics::plot(x$hclust,main="Hierarchical clustering of SHAP signatures",xlab="",sub=""); invisible(x) }
