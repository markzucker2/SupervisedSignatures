validate_required_columns <- function(x, required, object_name = "object") {
    if (!is.data.frame(x)) stop(object_name, " must be a data frame.")
    missing <- setdiff(required, names(x))
    if (length(missing)) stop(object_name, " is missing required column(s): ", paste(missing, collapse=", "))
    invisible(TRUE)
}

as_binary_outcome <- function(x, positive, negative) {
    factor(x[x %in% c(positive, negative)], levels=c(negative, positive))
}

encode_predictors <- function(x) {
    x <- as.data.frame(x, check.names=FALSE)
    for (nm in names(x)) if (is.character(x[[nm]])) x[[nm]] <- factor(x[[nm]])
    mm <- stats::model.matrix(~ . - 1, data=x)
    storage.mode(mm) <- "double"
    mm
}

safe_auc <- function(y, p) {
    y_num <- as.integer(y == levels(y)[2L])
    if (length(unique(y_num)) < 2L) return(NA_real_)
    as.numeric(pROC::auc(pROC::roc(y_num, p, quiet=TRUE)))
}

safe_auprc <- function(y, p) {
    y_num <- as.integer(y == levels(y)[2L])
    if (length(unique(y_num)) < 2L) return(NA_real_)
    as.numeric(PRROC::pr.curve(scores.class0=p[y_num==1], scores.class1=p[y_num==0], curve=FALSE)$auc.integral)
}

reshape_to_signature_matrix <- function(dat) {
    genes <- unique(dat$gene); features <- unique(dat$feature)
    mat <- matrix(0, nrow=length(genes), ncol=length(features), dimnames=list(genes,features))
    for (i in seq_len(nrow(dat))) mat[as.character(dat$gene[i]), as.character(dat$feature[i])] <- dat$mean_abs_shap[i]
    mat
}
