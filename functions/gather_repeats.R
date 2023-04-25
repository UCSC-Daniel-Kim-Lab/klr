gather_repeats <- function(
    .data, overwrite, genome, id_col="gene", overwrite_name="name",
    regions=c("Promoter_Region", "Junction", "3UTR"), overwrite_filter=NULL
)
{
    ## Required libraries
    require("tidyr")
    require("GenomicRanges")
    require("Biostrings")

    ## Name the regions
    regions = setNames(regions, regions)

    ## Filter overwrite if requested
    if (class(overwrite_filter)!="NULL")
    {
        if (class(overwrite_filter)=="expression")
        {
            overwrite = overwrite %>%
                filter(overwrite_filter)

        } else
        {
           stop("Stopping: bad filter") # Must use an expression for filtering
        }
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
                        filter(
                            eval(subsetter),
                            classification=reg
                        ) %>%
                        mutate(
                            repStop = repStart+(-1*repLength), repStrand
                        ) %>%
                        unite(c(overwrite_id, repName), col="name", sep='_') %>%
                        unite(c(name, repStart), col="instance", sep='_', remove=F) %>%
                        select(
                            chr=chrom,
                            start=repStart,
                            stop=repStop,
                            strand=repStrand,
                            name=instance,
                            gene_strand=genoStrand
                        ) %>%
                        filter(
                            !str_ends(chr, "_alt")
                        ) %>%
                        distinct()

                    fasta = Biostrings::getSeq(
                        genome, GenomicRanges::GRanges(bed)
                    )
                    names(fasta) = bed$name

                    return(list("bed"=bed, "fasta"=fasta))
                }
            )
        )
    } else
    { # Data is grouped - make seqs for each group
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
                        filter(
                            eval(subsetter),
                            classification==reg
                        ) %>%
                        mutate(
                            repStop = repStart+(-1*repLength), repStrand
                        ) %>%
                        unite(c(o_id, repName), col="name", sep='_') %>%
                        unite(c(name, repStart), col="instance", sep='_', remove=F) %>%
                        select(
                            chr=chrom,
                            start=repStart,
                            stop=repStop,
                            strand=repStrand,
                            name=instance,
                            gene_strand=genoStrand
                        ) %>%
                        filter(
                            !str_ends(chr, "_alt")
                        ) %>%
                        distinct()

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
            select(!c(".rows")) %>%
            unnest_longer(res) %>%
            rename(c(res_id="region", res="genome")) %>%
            unnest_longer(genome) %>%
            pivot_wider(names_from=genome_id, values_from=genome)

        return(groups)
    }


}