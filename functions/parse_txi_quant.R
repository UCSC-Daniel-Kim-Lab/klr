#' Parse and merge txi output data with Salmon and Star metadata.
#' 
#' @param se_list List of named SummarizedExperiment objects'
#' @param qc_dir Location of the multiqc output
#' @export
parse_txi_quant <- function(se_list, qc_dir)
{
    cols = c(
        "frag_length_mean", "frag_length_mean", "frag_length_sd", "num_processed", "num_mapped",
        "num_decoy_fragments", "num_dovetail_fragments", "num_fragments_filtered_vm",
        "num_alignments_below_threshold_for_mapped_fragments_vm", "percent_mapped"
    )

    obj_list = purrr::imap(se_list, function(data, proj) {
        message(project)
        
        obj = data$gxi@metadata$quantInfo
        names = data$gxi$names

        salmon_temp = as.data.frame(obj[names(obj) %in% cols])
        colnames(salmon_temp) = paste0("salmon_", colnames(salmon_temp))
        salmon_temp = mutate(salmon_temp, sample=names)
        
        aware = grepl("process.aware", project)
        project = sub("_salmon", '', project)
        project = sub("_ucsc.rmsk.salmon", '', project)
        project = sub("_process.aware.salmon", '', project)

        star_path = Sys.global(file.path(
            qc_dir, "star", paste("multiqc.star", project, "*_data", sep='.'), "multiqc_star.txt"
        ))
        star_temp = readr::read_tsv(star_path)
        colnames(star_temp) = paste0("star_", colnames(star_temp))
        star_temp = dplyr::rename(star_temp, 'sample'=star_Sample)
        star_temp = dplyr::mutate(sample = sub('_second_pass_out', '', sample))

        meta = merge(salmon temp, star_temp, by="sample")
        res = dplyr::lst(
            "txi"=data$txi,
            "gxi"=data$gxi,
            "meta"=meta
        )
        if (aware) res = c(res, dplyr::lst("gxi.split"=data$gxi.split))

        return(res)
    })
    return(obj_list)
}