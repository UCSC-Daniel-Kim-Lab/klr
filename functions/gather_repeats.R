#' Collect the DNA sequences of all repeats in a list of (grouped) genes
#'
#' @param .data A data frame, with a column of gene names or ensg names
#' @param overwrite A dataframe containing the repeat information for all genes
#' @param genome A genome object compatible with Biostrings::getSeq
#' @param id_col Name of the column in your dataframe holding the gene names
#' @param overwrite_name Name of teh column in your dataframe holding the gene names
#' @param overwrite_filter An expression to subset the overwite table to the desired set of repeats
#' @param regions The set of overwrite repeat classifications to include
#' @export
#' @example genes %>% group_by(expression) %>% gather_repeats(overwrite, BSgenome.hsapiens.hg38)
gather_repeats <- function(
    .data, overwrite, genome, id_col="gene", overwrite_name="name",
    regions=c("Promoter_Region", "Junction", "3UTR"), overwrite_filter=NULL
)
{
    ## Filter overwrite if requested
    if (class(overwrite_filter)!="NULL")
    {
        if (class(overwrite_filter)=="expression")
        {
            overwrite = overwrite %>%
                dplyr::filter(overwrite_filter)

        } else
        {
           stop("Stopping: bad filter") # Must use an expression for filtering
        }
    }

    ## Get the regions
    if (is.null(regions))
    {
        regions = unqiue(overwrite$classification)
        regions = stats::setNames(regions, regions)
    } else
    {
        regions = stats::setNames(regions, regions)
    }

    ## Generate Sequences
    if (is.null(attr(.data, "groups")))
    { # Data is not grouped - make all seqs
        # Get gene list
        genes = data[[id_col]]

        # Generate and return seqs
        return(
            imap(
                .x=regions, .f=function(reg, reg_name)
                {
                    # Subset overwrite
                    subsetter = expression(
                        overwrite[[overwrite_name]] %in% genes
                    )

                    bed = overwrite %>%
                        dplyr::filter(
                            eval(subsetter),
                            classification=reg
                        ) %>%
                        dplyr::mutate(
                            repStop = repStart+repLength, repStrand
                        ) %>%
                        tidyr::unite(c(overwrite_id, repName), col="name", sep='_') %>%
                        tidyr::unite(c(name, repStart), col="instance", sep='_', remove=F) %>%
                        dplyr::select(
                            chr=chrom,
                            start=repStart,
                            stop=repStop,
                            strand=repStrand,
                            name=instance,
                            gene_strand=genoStrand
                        ) %>%
                        dplyr::filter(
                            !str_ends(chr, "_alt")
                        ) %>%
                        dplyr::distinct()

                    fasta = Biostrings::getSeq(
                        genome, GenomicRanges::GRanges(bed)
                    )
                    names(fasta) = bed$name

                    return(list("bed"=bed, "fasta"=fasta))
                }
            )
        )
    }
    # Data is grouped - make seqs for each group
    # Gather grouped data and grouping vars
    groups = attr(.data, "groups")
    cols = colnames(groups)[colnames(groups) != ".rows"]

    # Make a 'vectorized' version of the above
    seqs_by_row = function(rows, o, g, id, o_id, regs)
    {
        # Get list of genes
        genes = .data[rows, ][[id_col]]

        # Gather repeats
        seqs = imap(
            .x=regions, .f=function(rg, reg_name)
            {
                # Make overwrite filter
                subsetter = expression(
                    o[[o_id]] %in% genes
                )

                bed = o %>%
                    dplyr::filter(
                        eval(subsetter),
                        classification==reg
                    ) %>%
                    dplyr::mutate(
                        repStop = repStart+(-1*repLength), repStrand
                    ) %>%
                    tidyr::unite(c(o_id, repName), col="name", sep='_') %>%
                    tidyr::unite(c(name, repStart), col="instance", sep='_', remove=F) %>%
                    dplyr::select(
                        chr=chrom,
                        start=repStart,
                        stop=repStop,
                        strand=repStrand,
                        name=instance,
                        gene_strand=genoStrand
                    ) %>%
                    dplyr::filter(
                        !str_ends(chr, "_alt")
                    ) %>%
                    dplyr::distinct()

                fasta = Biostrings::getSeq(
                    g, GenomicRanges::GRanges(bed)
                )
                names(fasta) = bed$name

                return(list("bed"=bed, "fasta"=fasta))
            }
        )
        return(seqs)
    }

    ## Add DNAStringSet to the grouped data
    groups$res = groups %>%
        apply(MARGIN=1, FUN=function(row)
            {seqs_by_row(row$.row, overwrite, genome, id_col, overwrite_name, regions)}
        )
    
    ## Unnest results
    groups = groups %>%
        dplyr::select(!c(".rows")) %>%
        tidyr::unnest_longer(res) %>%
        dplyr::rename(c(res_id="region", res="genome")) %>%
        tidyr::unnest_longer(genome) %>%
        tidyr::pivot_wider(names_from=genome_id, values_from=genome)

    return(groups)
}