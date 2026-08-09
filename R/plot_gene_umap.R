#' Plot genes in SHAP signature space using UMAP
#' @param results A `SignatureResults` object.
#' @param genes Optional genes to include.
#' @param comparison Optional comparison ID.
#' @param n_neighbors UMAP neighborhood parameter.
#' @param min_dist UMAP minimum distance parameter.
#' @param metric UMAP metric.
#' @param seed Random seed.
#' @return UMAP coordinates invisibly; prints the plot.
#' @export
plot_gene_umap <- function(results,genes=NULL,comparison=NULL,n_neighbors=15L,min_dist=0.1,metric="cosine",seed=1L) {
    if (!inherits(results,"SignatureResults")) stop("results must be a SignatureResults object.")
    if (!requireNamespace("uwot",quietly=TRUE)) stop("Package 'uwot' is required for plot_gene_umap().")
    dat <- results$shap_summary; if (!is.null(genes)) dat <- dat[dat$gene %in% genes,,drop=FALSE]; if (!is.null(comparison)) dat <- dat[dat$comparison_id==comparison,,drop=FALSE]
    if (!nrow(dat)) stop("No SHAP summary data available.")
    gene_counts <- tapply(dat$comparison_id,dat$gene,function(x) length(unique(x))); if (any(gene_counts>1L)) stop("Multiple comparisons exist for at least one gene. Specify comparison.")
    mat <- reshape_to_signature_matrix(dat); if (nrow(mat)<3L) stop("At least 3 genes are required for UMAP.")
    set.seed(seed); embedding <- uwot::umap(mat,n_neighbors=min(n_neighbors,nrow(mat)-1L),min_dist=min_dist,metric=metric,verbose=FALSE)
    coords <- data.frame(gene=rownames(mat),UMAP1=embedding[,1],UMAP2=embedding[,2],stringsAsFactors=FALSE)
    p <- ggplot2::ggplot(coords,ggplot2::aes(UMAP1,UMAP2,label=gene)) + ggplot2::geom_point() + ggplot2::geom_text(vjust=-0.7,check_overlap=TRUE) +
        ggplot2::labs(title="Genes in SHAP signature space",x="UMAP 1",y="UMAP 2") + ggplot2::theme_bw()
    print(p); invisible(coords)
}
