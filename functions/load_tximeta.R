#' Read and merge the transcript and gene level SE objects.
#' 
#' @param se_dir The location of the SE objects.
#' @param ref The 
#' 
load_tximeta <- function(se_dir, ref)
{
    out_path = file.path(se_dir, paste0(ref, "_quant"))
    obj_list = list.files(se_dir)
    temp = sub("_h5_se", '', obj_list)
    temp = sub("_gene", '', temp)
    temp = sub("_split", '', temp)
    sample_list = unique(temp)
    names(sample_list) = sample_list
    print(sample_list)

    sample_list %>%
        purrr::imap(function(name, sample)
        {
            txi = HDF5Array::loadHDF5SummarizedExperiment(
                dir=file.path(se_dir, paste0(sample, "_h5_se"))
            )
            gxi = HDF5Array::loadHDF5SummarizedExperiment(
                dir=file.path(se_dir, paste0(sample, "_gene_h5_se"))
            )

            if (ref=="process.aware.salmon")
            {
                gxi_split = HDF5Array::loadHDF5SummarizedExperiment(
                    dir=file.path(se_dir, paste0(sample, "_gene_split_h5_se"))
                )
                return(list("txi"=txi, "gxi"=gxi, "gxi.split"=gxi_split))
            } else
            {
               return(list("txi"=txi, "gxi"=gxi))
            }
        }) 
}