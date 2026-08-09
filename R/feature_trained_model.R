### old version; deprecated
library(data.table)
library(caret)
library(xgboost)
library(pROC)
library(PRROC)
library(skitools)
library(plyr)
library(SHAPforxgboost, lib="/gpfs/home/zuckem04/projects/CIN/scripts/R_libPath")
set.seed(123)

download.file("https://cran.r-project.org/src/contrib/Archive/MLwrap/MLwrap_0.1.0.tar.gz",
              "~/git/MLwrap_0.1.0.tar.gz")
install.packages("~/git/MLwrap_0.1.0.tar.gz", type='source', repos=NULL,
                 lib="/gpfs/home/zuckem04/projects/CIN/scripts/R_libPath")

#Am I local or am I on the cluster:
cd <- getwd()
if(grepl('gpfs',cd)){data_path <- '~/projects/CIN/files/data'}else{data_path <- '/Users/zuckem04/Documents/CIN/models/data'}


### testing:
#testing:
features_hmf <- readRDS(paste0(data_path, '/hmf/merged_features_hmf.rds'))
features_tcga <- readRDS(paste0(data_path, '/tcga/merged_features_tcga.rds'))
features_nygc <- readRDS(paste0(data_path, '/nygc/merged_features_nygc.rds'))
features_all <- rbind(features_hmf, features_nygc)
colnames(features_all)[1] <- 'sample_id'
gt <- readRDS(paste0(data_path, '/merged/merged_gt.rds'))
features_all <- rbind(features_hmf, features_nygc, features_tcga)
colnames(features_all)[1] <- 'sample_id'
#converting to numeric:
for(i in 2:ncol(features_all)){features_all[,i] <- as.numeric(features_all[,i])}

#signatures:
merged_signatures <- readRDS(paste0(data_path, '/merged/merged_signatures_all.rds'))

#test run on MAD1L1, BUB1B, BRCA1, BRCA2
g <- c('MAD1L1','MAD2L1','BUB1B','BRCA1','BRCA2','TP53')
gene_states_all <- gt[gene %in% g, c('pair','gt','gene')]
colnames(gene_states_all) <- c('sample_id', 'state','gene')
gene_states_all$state[gene_states_all$state %in% c('homdel','mutloh')] <- "biallelic"

#Load pairs_merged to get covariates:
pairs_merged <- readRDS(paste0(data_path, '/merged/pairs_merged.rds'))
colnames(pairs_merged)[c(2,5,7)] <- c('tumor_type_old','metastatic_status','tumor_type')
colnames(pairs_merged)[1] <- 'sample_id'
samples_all <- pairs_merged[,c(1, 3:7)]
samples_all$metastatic_status <- mapvalues(samples_all$metastatic_status, 
                                           from=c('TRUE','FALSE'), to=c('metastatic','primary'))
samples_all$ploidy <- as.numeric(samples_all$ploidy)
samples_all$purity <- as.numeric(samples_all$purity)
rm(gt, g, features_hmf, features_tcga, features_nygc)

#saving if files don't already exist
sample_path <- paste0(data_path, '/samples.rds')
feature_path <- paste0(data_path, '/features.rds')
genestate_path <- paste0(data_path, '/gene_states.rds')
if(!file.exists(sample_path)){saveRDS(samples_all, file=sample_path)}
if(!file.exists(feature_path)){saveRDS(features_all, file=feature_path)}
if(!file.exists(genestate_path)){saveRDS(gene_states_all, file=genestate_path)}

#for now, exlcuding TCGA/hg38 patients:
#dat <- dat[dataset != 'TCGA8.5k']

#First, let's check b1/b2, other sanity checks
#b1 <- gt[gene == 'BRCA1']
#b1$B1 <- merged_signatures$B1[match(b1$pair, merged_signatures$pair)]
#ppdf(print(ggplot(data=b1, aes(x=gt, y=B1, fill=gt)) + geom_violin()))
#b2 <- gt[gene == 'BRCA2']
#b2$B2 <- merged_signatures$B2[match(b2$pair, merged_signatures$pair)]
#ppdf(print(ggplot(data=b2, aes(x=gt, y=B2, fill=gt)) + geom_violin()))
#palb2 <- gt[gene == 'PALB2']
#palb2$hrdetect <- merged_signatures$hrdetect[match(palb2$pair, merged_signatures$pair)]
#ppdf(print(ggplot(data=palb2, aes(x=gt, y=hrdetect, fill=gt)) + geom_violin()))
#cdkn2a <- gt[gene == 'CDKN2A']
#table(cdkn2a$gt[cdkn2a$tumor_type == ''])

#Load data:
features_all <- readRDS(sample_path)
samples <- readRDS(feature_path)
gene_states_all <- readRDS(genestate_path)

###################################################
# Pipeline
###################################################

STATE_LEVELS <- c(
  "wt",
  "hetloss",
  "biallelic",
  "gain",
  "amp"
)

CLASSIFIERS <- list(
  list(name="hetloss_vs_wt",
       positive=c("hetloss"),
       negative=c("wt")),
  list(name="biallelic_vs_wt",
       positive=c("biallelic"),
       negative=c("wt")),
  list(name="gain_vs_wt",
       positive=c("gain"),
       negative=c("wt")),
  list(name="amp_vs_wt",
       positive=c("amp"),
       negative=c("wt"))
)

# Covariates
COVARIATE_COLUMNS <- c(
  "metastatic_status",
  "dataset",
  "ploidy",
  "purity"
)
#to add: TMB, SV_TMB, MSI_status; tumor_type?

FEATURE_COLUMNS <- colnames(features_all[,-1])

validate_inputs <- function(samples,features,gene_states){
  stopifnot("sample_id"%in%names(samples),"sample_id"%in%names(features))
  stopifnot(all(c("sample_id","gene","state")%in%names(gene_states)))
}

build_dataset<-function(gene_name,positive_state,negative_state,samples,features,gene_states){
  gs <- gene_states[gene==gene_name & state %in% c(positive_state,negative_state)]
  gs[,y:=ifelse(state==positive_state,1L,0L)]
  dat <- merge(samples,features,by="sample_id")
  dat <- merge(dat,gs[,.(sample_id,y)],by="sample_id")
  if("TMB"%in%names(dat)) dat[,TMB:=log10(TMB+1)]
  if("SV_TMB"%in%names(dat)) dat[,SV_TMB:=log10(SV_TMB+1)]
  dat
}

fit_classifier <- function(dat){
  feature_cols <- FEATURE_COLUMNS[FEATURE_COLUMNS %in% names(dat)]
  predictor_cols <- c(COVARIATE_COLUMNS[COVARIATE_COLUMNS%in%names(dat)],feature_cols)
  idx <- createDataPartition(dat$y,p=.8,list=FALSE)
  train <- dat[idx]; test <- dat[-idx]
  
  #ensuring levels are the same for training and test model matrices
  all_levels <- union(unique(train$tumor_type), unique(test$tumor_type))
  train[, tumor_type := factor(tumor_type, levels = all_levels)]
  test[, tumor_type := factor(tumor_type, levels = all_levels)]
  
  #Making model matrices
  xtrain <- model.matrix(~.-1,data=train[,..predictor_cols]); xtest <- model.matrix(~.-1,data=test[,..predictor_cols])
  dtrain <- xgb.DMatrix(xtrain,label=train$y); dtest <- xgb.DMatrix(xtest,label=test$y)
  model <- xgb.train(params=list(objective="binary:logistic",eval_metric="auc",
                                 eta=.05,max_depth=6,subsample=.8,colsample_bytree=.8,min_child_weight=5),
                     data=dtrain,nrounds=1000,watchlist=list(train=dtrain,test=dtest),
                     early_stopping_rounds=50,verbose=1)
  pred <- predict(model,dtest)
  rocobj <- roc(test$y,pred,quiet=TRUE)
  pr <- pr.curve(scores.class0=pred[test$y==1],scores.class1=pred[test$y==0])
  
  #shap values:
  shap <- shap.values(
    xgb_model=model,
    X_train=xtrain
  )
  
  shap_long <- shap.prep(
    shap_contrib=shap$shap_score,
    X_train=xtrain
  )
  
  list(model=model,AUROC=as.numeric(auc(rocobj)),AUPRC=pr$auc.integral,
       importance=xgb.importance(colnames(xtrain),model),
       shap=shap, shap_long=shap_long)
}

run_pipeline <- function(samples,features,gene_states,classifiers=CLASSIFIERS){
  validate_inputs(samples,features,gene_states)
  genes <- sort(unique(gene_states$gene));fits<-list();summary<-data.table()
  
  for(g in genes){
    for(cl in classifiers){
      dat <- build_dataset(g,cl$positive,cl$negative,samples,features,gene_states)
      dat <- na.omit(dat)
      if(nrow(dat)==0||length(unique(dat$y))<2) next
      fit <- fit_classifier(dat)
      fits[[paste(g,cl$name,sep="_")]] <- fit
      summary <- rbind(summary,data.table(gene=g,comparison=cl$name,AUROC=fit$AUROC,AUPRC=fit$AUPRC,n_positive=sum(dat$y==1),n_negative=sum(dat$y==0)))
    }}
  list(summary=summary,models=fits)
}


## Example:
#excluding rare cancer types:
#typetbl <- sort(table(samples$tumor_type), decreasing=T)
#samples <- samples_all[samples_all$tumor_type %in% names(typetbl[typetbl >= 100]),]
#features <- features_all[match(features$sample_id,samples_all$sample_id),]
#gene_states <- gene_states_all[match(gene_states_all$sample_id,samples_all$sample_id),]
samples <- samples_all
features <- features_all
gene_states <- gene_states_all

#Restricting features to just SV, CN, and event features for now:
features <- features[,grepl('sample_id', colnames(features)) | 
                       !grepl('Del|Ins|seg|LOOSE', colnames(features))]

result <- run_pipeline(samples,features,gene_states)
result$summary
saveRDS(result, file='/gpfs/home/zuckem04/projects/CIN/flow_supervised_model/results/results_7_17_26.rds')
plot_shap(result,gene="BRCA2",state="LOSS")

## Run on validation set (TCGA data)



### Need good neg. controls: 
#




