test_that("SignatureDataset validates basic inputs", {
    samples <- data.frame(sample_id=paste0("S",1:3))
    features <- data.frame(sample_id=paste0("S",1:3), f1=1:3, f2=3:1)
    covariates <- data.frame(sample_id=paste0("S",1:3), tumor_type=c("A","B","A"))
    gene_states <- data.frame(sample_id=paste0("S",1:3), gene="BRCA2", state=c("WT","LOSS","MUT"))
    ds <- SignatureDataset(samples,features,covariates,gene_states)
    expect_s3_class(ds,"SignatureDataset")
    expect_equal(ncol(ds$features)-1,2)
})

test_that("comparison generator creates altered-vs-WT comparisons", {
    gs <- data.frame(sample_id=paste0("S",1:6),gene=rep(c("BRCA1","BRCA2"),each=3),state=rep(c("WT","LOSS","MUT"),2))
    cmp <- make_gene_state_comparisons(gs)
    expect_equal(nrow(cmp),4)
    expect_true(all(cmp$negative=="WT"))
})

test_that("build_dataset keeps feature columns separate from covariates", {
    samples <- data.frame(sample_id=paste0("S",1:4))
    features <- data.frame(sample_id=paste0("S",1:4),SV1=c(1,2,3,4),CN1=c(4,3,2,1))
    covariates <- data.frame(sample_id=paste0("S",1:4),tumor_type=c("A","A","B","B"))
    gene_states <- data.frame(sample_id=paste0("S",1:4),gene="BRCA2",state=c("WT","WT","LOSS","LOSS"))
    ds <- SignatureDataset(samples,features,covariates,gene_states)
    cmp <- make_gene_state_comparisons(gene_states)
    ad <- build_dataset(ds,cmp[1,])
    expect_equal(ad$feature_columns,c("SV1","CN1"))
    expect_equal(ad$covariate_columns,"tumor_type")
    expect_equal(ncol(ad$X),3)
})
