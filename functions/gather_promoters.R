gather_promoters <- function(genes, genome, gtf_file, id_col="gene")
{
    if (is.null(gtf_path) | is.null(genome))
    {
        stop("Must provide a BSGenome and GTF file path")
    }

    gtf = rtracklayer::readGFFAsGRanges(filepath=gtf_file)

    if (is.null(genes))
    {
        gtf %>%
            as.data.frame() %>%
            filter(type="gene") %>%
            GenomicRanges::GRanges() %>%
            GenomicRanges::promoters()
        seqnames = gtf$gene_id
        gtf = BioStrings::getSeq(genome, gtf)
        names(gtf) = seqnames
        return(gtf)
    }
    else if ("vector" %in% class(genes) | "list" %in% class(genes))
    {
        gtf %>%
            as.data.frame() %>%
            filter(
                type=="gene",
                gene_id %in% genes
            ) %>%
            GenomicRanges::GRanges() %>%
            GenomicRanges::promoters()
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

            gtf %>%
                as.data.frame() %>%
                filter(
                    type=="gene",
                    gene_id %in% gene_list
                ) %>%
                GenomicRanges::GRanges() %>%
                GenomicRanges::promoters()
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
                bed %>%
                    as.data.frame() %>%
                    filter(type=="gene") %>%
                    GenomicRanges::GRanges() %>%
                    GenomicRanges::promoters()
                seqnames = gtf$gene_id
                bed = BioStrings::getSeq(genome, bed)
                names(bed) = seqnames
                return(gtf)
            }

            groups$promoters <- unlist(apply(
                groups, MARGIN=1, FUN=promoters_by_row,
                data=genes, id_col=id_col, bed=gtf, genome=genome
            ))
            
            groups %>%
                dplyr::select(!c(".rows")) %>%
                return()
        }
    }
}