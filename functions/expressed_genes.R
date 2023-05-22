#' Get a list of all expressed genes from Salmon Quantification
#' 
#' @param h5_path
#' @param tx_to_gene
#' @param id_col
#' @param gene_col
#' @param raw
#' @param count_threshohld
#' @param count_freq
#' @param raw
#' @param biotyper
expressed_genes <- function(
    h5_path, tx_to_gene, id_col="gene", gene_col="V2"
    raw=F, count_threshhold=5, count_freq=0.75, biotyper=NULL,
)
{
    biotypes = dplyr::select(
        dplyr::ungroup(
            dplyr::distinct(
                dplyr::filter(
                    dplyr::mutate(
                        dplyr::group_by(
                            dplyr::select(
                                tx_to_gene,
                                gene_col, biotype
                            ),
                            gene_col
                        ),
                        biotype=case_when(
                            all(biotype=="lncRNA")~"lncRNA",
                            any(biotype=="protein_coding")~"protein_coding",
                            .default="remove"
                        )
                    ),
                    biotype!="remove"
                )
            )
        ),
        biotype, gene_col
    )
    se = HDF5Array::loadHDF5SummarizedExperiment(h5_path)
    counts = {}
    if (raw)
    {
        counts = SummarizedExperiment::assays(se)$counts
    } else
    {
        counts = SummarizedExperiment::assays(se)$abundance
    }

    genes = rownames(counts[(rowSums(counts >= count_threshhold)) >= ncol(counts)*count_freq, ])
    genes = merge(
            x=data.frame(gene=genes), y=biotypes,
            by.x="gene", by.y=gene_col
        )
    return(genes)
}