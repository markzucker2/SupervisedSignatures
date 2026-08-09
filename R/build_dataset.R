#' Build an analysis dataset for one gene-state comparison
#' @param dataset A `SignatureDataset`.
#' @param comparison One row of a comparison data frame containing `gene`, `positive`, and `negative`.
#' @return An object of class `SignatureAnalysisDataset`.
#' @export
build_dataset <- function(dataset, comparison) {
    if (!inherits(dataset,"SignatureDataset")) stop("dataset must be a SignatureDataset.")
    if (!all(c("gene","positive","negative") %in% names(comparison))) stop("comparison must contain gene, positive, and negative.")
    gene <- as.character(comparison$gene[[1L]]); positive <- as.character(comparison$positive[[1L]]); negative <- as.character(comparison$negative[[1L]])
    comparison_id <- if ("id" %in% names(comparison)) as.character(comparison$id[[1L]]) else paste(gene,positive,"vs",negative,sep="_")
    gs <- dataset$gene_states[dataset$gene_states$gene==gene & dataset$gene_states$state %in% c(positive,negative), c("sample_id","gene","state"), drop=FALSE]
    if (!nrow(gs)) stop("No gene-state observations found for ",gene,".")
    dat <- merge(dataset$samples,dataset$features,by="sample_id",all=FALSE)
    if (!is.null(dataset$covariates)) dat <- merge(dat,dataset$covariates,by="sample_id",all=FALSE)
    dat <- merge(dat,gs,by="sample_id",all=FALSE)
    dat <- dat[dat$state %in% c(positive,negative),,drop=FALSE]
    feature_cols <- setdiff(names(dataset$features),"sample_id")
    covariate_cols <- if (is.null(dataset$covariates)) character() else setdiff(names(dataset$covariates),"sample_id")
    x_raw <- dat[,c(feature_cols,covariate_cols),drop=FALSE]
    y <- as_binary_outcome(dat$state,positive,negative)
    if (length(unique(y)) != 2L) stop("The selected comparison requires both outcome classes.")
    out <- list(data=dat,X_raw=x_raw,X=encode_predictors(x_raw),y=y,sample_id=dat$sample_id,
                gene=gene,positive=positive,negative=negative,comparison_id=comparison_id,
                feature_columns=feature_cols,covariate_columns=covariate_cols)
    class(out) <- "SignatureAnalysisDataset"; out
}
