#' Get a list of all expressed genes from Salmon Quantification
#' 
#' @param h5_path
expressed_genes <- function(
    h5_path, tx_to_gene, id_col="gene", gene_col="V2"
    raw=F, count_threshold=5, count_freq=0.75, biotyper=NULL,
)
{
    biotypes = tx_to_gene %>%
        dplyr::select(gene_col, biotype) %>%
        dplyr::group_by(gene_col) %>%
        dplyr::mutate(
            biotype = dplyr::case_when(
                all(biotype=="lncRNA") ~ "lncRNA",
                any(biotype=="protein_coding") ~ "protein_coding",
                .default="remove"
            )
        ) %>%
        dplyr::filter(
            biotype != "remove"
        ) %>%  
        dplyr::distinct() %>%
        dplyr::ungroup() %>%
        dplyr::select(biotype, gene_col)

    se = HDF5Array::loadHDF5SummarizedExperiment(h5_path)
    counts = {}
    if (raw)
    {
        counts = SummarizedExperiment::assays(se)$counts
    } else
    {
        counts = SummarizedExperiment::assays(se)$abundance
    }

    genes = rownames(counts[(rowSums(counts >= count_threshold)) >= ncol(counts)*count_freq, ]) %>%
        merge(
            x=data.frame(gene=.), y=biotypes,
            by.x="gene", by.y=""
        )
    return(genes)
}