# SupervisedSignatures

Supervised genomic signature discovery using machine learning.

## Initial development version

This version supports `SignatureDataset()`, automatic gene × altered-state-vs-reference comparisons, binary XGBoost classifiers, train/test AUROC and AUPRC, native XGBoost SHAP values, SHAP summary barplots, UMAP of genes in SHAP signature space, and hierarchical clustering.

```r
library(SupervisedSignatures)

dataset <- SignatureDataset(
    samples = samples,
    features = signature_features,
    covariates = covariates,
    gene_states = gene_states,
    feature_type = "COSMIC",
    feature_set = "SV+CN"
)

results <- discover_signatures(dataset, reference_state="WT", compute_shap=TRUE)
plot_shap(results, gene="BRCA2")
coords <- plot_gene_umap(results)
clusters <- cluster_signatures(results, k=5)
plot(clusters)
```

## Input conventions

`features` must contain `sample_id`. Every other column is treated as a model feature. `covariates`, if supplied, must contain `sample_id`; every other column is treated as a covariate and categorical variables are encoded automatically. `gene_states` must contain `sample_id`, `gene`, and `state`.

This is an early development version; APIs will evolve.
