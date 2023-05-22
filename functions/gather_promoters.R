#' Collect the promoter sequences from a set of genes using a BSgenome and an annotation
#'
#' @param genes The gene set as Ensemble Gene IDs. Can be a vector or data frame. Returns one list
#'  of promoter sequences unless the data frame is grouped, in which case it returns one list for
#'  each group.
#' @param genome The BSGenome object from which to pull the DNA sequences
#' @param gtf_file The genome annotation which contains the promoter regions.
#' @param id_col The column containing the Gene IDs, if a data frame is used.
#' @export
gather_promoters <- function(genes, genome, gtf_file, id_col="gene")
{
    if (is.null(gtf_path) | is.null(genome))
    {
        stop("Must provide a BSGenome and GTF file path")
    }

    gtf = rtracklayer::readGFFAsGRanges(filepath=gtf_file)

    if (is.null(genes))
    {
        gtf = GenomicRanges::promoters(
            GenomicRanges::GRanges(
                filter(
                    as.data.frame(gtf),
                    type=="gene"
                )
            )
        )
        seqnames = gtf$gene_id
        gtf = BioStrings::getSeq(genome, gtf)
        names(gtf) = seqnames
        return(gtf)
    }
    else if ("vector" %in% class(genes) | "list" %in% class(genes))
    {
        gtf = GenomicRanges::promoters(
            GenomicRanges::GRanges(
                filter(
                    as.data.frame(gtf),
                    type=="gene",
                    gene_id=="gene"
                )
            )
        )
        seqnames = gtf$gene_id
        gtf = BioStrings::getSeq(genome, gtf) %>%
        names(gtf) = seqnames
        return(gtf)
    }
    else if ("data.frame" %in% class(genes))
    {
        if (is.null(attr(genes, 'groups')))
        {
            gene_list = genes[[id_col]]
            gtf = GenomicRanges::promoters(
                GenomicRanges::GRanges(
                    filter(
                        as.data.frame(gtf),
                        type=="gene",
                        gene_id %in% gene_list
                    )
                )
            )
            seqnames = gtf$gene_id
            gtf = BioStrings::getSeq(genome, gtf)
            names(gtf) = seqnames
            return(gtf)
        }
        else
        {
            groups = attr(genes, 'groups')
            cols = colnames(groups)[colnames(groups)!=".row"]

            promoters_by_row = function(group, data, id_col, bed, genome)
            {
                gene_list = data[[id_col]][group$.rows]
                bed = GenomicRanges::promoters(
                    GenomicRanges:: GRanges(
                        filter(
                            as.data.frame(bed),
                            type=="gene"
                        )
                    )
                )
                seqnames = gtf$gene_id
                bed = BioStrings::getSeq(genome, bed)
                names(bed) = seqnames
                return(gtf)
            }

            groups$promoters <- unlist(apply(
                groups, MARGIN=1, FUN=promoters_by_row,
                data=genes, id_col=id_col, bed=gtf, genome=genome
            ))
            
            return(dplyr::select(groups, !c(".rows")))
        }
    }
}