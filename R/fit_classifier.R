#' Fit a binary XGBoost classifier
#' @param analysis_dataset A `SignatureAnalysisDataset`.
#' @param test_fraction Fraction assigned to the test set.
#' @param seed Random seed.
#' @param params List of XGBoost parameters.
#' @param nrounds Number of boosting rounds.
#' @param early_stopping_rounds Optional early stopping rounds.
#' @return A `SignatureFit` object.
#' @export
fit_classifier <- function(analysis_dataset,test_fraction=0.20,seed=1L,
    params=list(objective="binary:logistic",eval_metric="auc",max_depth=4,eta=0.05,subsample=0.8,colsample_bytree=0.8,min_child_weight=1),
    nrounds=500L,early_stopping_rounds=30L) {
    if (!inherits(analysis_dataset,"SignatureAnalysisDataset")) stop("analysis_dataset must be a SignatureAnalysisDataset.")
    n <- nrow(analysis_dataset$X); if (n < 10L) stop("At least 10 samples are recommended.")
    set.seed(seed); idx <- sample.int(n); n_test <- max(1L,floor(n*test_fraction)); test_idx <- idx[seq_len(n_test)]; train_idx <- idx[-seq_len(n_test)]
    X_train <- analysis_dataset$X[train_idx,,drop=FALSE]; X_test <- analysis_dataset$X[test_idx,,drop=FALSE]
    y_train <- as.integer(analysis_dataset$y[train_idx])-1L; y_test <- as.integer(analysis_dataset$y[test_idx])-1L
    dtrain <- xgboost::xgb.DMatrix(X_train,label=y_train); dtest <- xgboost::xgb.DMatrix(X_test,label=y_test)
    args <- list(data=dtrain,nrounds=nrounds,params=params,watchlist=list(train=dtrain,test=dtest),verbose=0)
    if (!is.null(early_stopping_rounds)) args$early_stopping_rounds <- early_stopping_rounds
    model <- do.call(xgboost::xgb.train,args)
    pred_train <- as.numeric(stats::predict(model,dtrain)); pred_test <- as.numeric(stats::predict(model,dtest))
    ytr <- factor(y_train,levels=c(0,1),labels=c(analysis_dataset$negative,analysis_dataset$positive)); yte <- factor(y_test,levels=c(0,1),labels=c(analysis_dataset$negative,analysis_dataset$positive))
    metrics <- data.frame(comparison_id=analysis_dataset$comparison_id,gene=analysis_dataset$gene,positive=analysis_dataset$positive,negative=analysis_dataset$negative,
        n_train=length(y_train),n_test=length(y_test),prevalence_train=mean(y_train),prevalence_test=mean(y_test),
        auc_train=safe_auc(ytr,pred_train),auc_test=safe_auc(yte,pred_test),auprc_train=safe_auprc(ytr,pred_train),auprc_test=safe_auprc(yte,pred_test),stringsAsFactors=FALSE)
    importance <- xgboost::xgb.importance(feature_names=colnames(X_train),model=model)
    out <- list(model=model,train_idx=train_idx,test_idx=test_idx,X_train=X_train,X_test=X_test,y_train=ytr,y_test=yte,
        predictions=list(train=pred_train,test=pred_test),metrics=metrics,importance=importance,comparison_id=analysis_dataset$comparison_id,
        gene=analysis_dataset$gene,positive=analysis_dataset$positive,negative=analysis_dataset$negative,sample_id_train=analysis_dataset$sample_id[train_idx],
        sample_id_test=analysis_dataset$sample_id[test_idx],feature_columns=analysis_dataset$feature_columns,covariate_columns=analysis_dataset$covariate_columns)
    class(out) <- "SignatureFit"; out
}
