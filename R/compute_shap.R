#' Compute native XGBoost SHAP values
#' @param fit A `SignatureFit`.
#' @return The input fit with `shap_raw` and `shap_summary` added.
#' @export
compute_shap <- function(fit) {
    if (!inherits(fit,"SignatureFit")) stop("fit must be a SignatureFit.")
    shap <- stats::predict(fit$model,xgboost::xgb.DMatrix(fit$X_train),predcontrib=TRUE); shap <- as.matrix(shap)
    feature_names <- colnames(fit$X_train)
    if (ncol(shap)==length(feature_names)+1L) { colnames(shap) <- c(feature_names,"BIAS"); shap_features <- shap[,feature_names,drop=FALSE] }
    else if (ncol(shap)==length(feature_names)) { colnames(shap) <- feature_names; shap_features <- shap }
    else stop("Unexpected number of SHAP columns returned by XGBoost.")
    mean_abs <- colMeans(abs(shap_features),na.rm=TRUE)
    summary <- data.frame(comparison_id=fit$comparison_id,gene=fit$gene,positive=fit$positive,negative=fit$negative,
                          feature=names(mean_abs),mean_abs_shap=as.numeric(mean_abs),stringsAsFactors=FALSE)
    summary <- summary[order(summary$mean_abs_shap,decreasing=TRUE),,drop=FALSE]
    fit$shap_raw <- shap; fit$shap_summary <- summary; fit
}
