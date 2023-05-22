#' Widen the motifs (and optionally FIMO sequences) into a much larger, readable data frame.
#' 
#' @param motifs The dataframe containing the grouped motif objects.
#' @param id_col The column in the data frame containing the grouped motifs
#' @param var_cols The columns containing important metadata about the grouped motifs.
#' @param fimo [OPTIONAL] The column containing the FIMO output results from the motifs.
#' @export
expand_motifs <- function(motifs, id_col="motif", var_cols, fimo)
{
    cols =  c(id_col)
    if (!is.null(fimo)) cols = c(cols, fimo)
    if (!is.null(var_cols)) cols = c(cols, var_cols)
    
    motifs = tidyr::unnest(
        dplyr::select(motifs, cols),
        cols=tidyr::matches(id_col)
    )
    
    if (!is.null(fimo))
    {
        motifs = tidyr::unnest_longer(
            motifs,
            cols=tidyr::matches(fimo)
        )
    }

    return(motifs)
}