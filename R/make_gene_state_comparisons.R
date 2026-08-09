#' Generate default gene-state comparisons
#'
#' @param gene_states Data frame with `gene` and `state` columns.
#' @param reference_state Reference state, usually `"WT"`.
#' @param states Optional vector of states to include. If `NULL`, all states other than the reference are used.
#' @return A data frame with `id`, `gene`, `positive`, and `negative`.
#' @export
make_gene_state_comparisons <- function(gene_states, reference_state="WT", states=NULL) {
    validate_required_columns(gene_states, c("gene","state"), "gene_states")
    genes <- unique(as.character(gene_states$gene)); out <- list(); k <- 0L
    for (g in genes) {
        available <- unique(as.character(gene_states$state[gene_states$gene==g]))
        selected <- if (is.null(states)) setdiff(available, reference_state) else intersect(states, available)
        for (s in selected) {
            k <- k + 1L
            out[[k]] <- data.frame(id=paste(g,s,"vs",reference_state,sep="_"), gene=g,
                                   positive=s, negative=reference_state, stringsAsFactors=FALSE)
        }
    }
    if (!length(out)) return(data.frame(id=character(),gene=character(),positive=character(),negative=character(),stringsAsFactors=FALSE))
    do.call(rbind,out)
}
