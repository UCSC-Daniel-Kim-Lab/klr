#' Find the enriched set of motifs in a set of DNA sequences using an existing database of motifs.
#' NOTE: This requires the installation of the MEME Suite _*AND*_ the memes R library.
#' 
#' @param .data The set of DNA strings. A single list, or a data frame object can be used.
#' @param .y A second list of DNA strings to be used as the control or background (not used with the
#'  data frame)
#' @param database The path to the meme formatted database file, or a list of paths to many meme
#'  formatted databases.
#' @param genes_col The col containing the DNA sequences to search for motifs. (Only used if .data
#'  is a data.frame)
#' @param background_col The column containing the DNA sequences for the background genes. (Only 
#'  used if .data is a data.frame)
#' @export
locate_motifs <- function(.data, .y, database, genes_col="obs", background_col="bg")
{
    res = NULL
    if (class(.data)=="data.frame")
    {
        res = unlist(apply(
            .data, MARGIN=1, FUN=function(row)
            {
                temp = try(memes::runAme(
                    row[[genes_col]], control=row[[background_col]],
                    database=database
                ))
                if (is.null(temp) | "try-error" %in% class(temp))
                {
                    return("No motifs found")
                }
                return(temp)
            }
        ))
    }
    else if (class(.data)=="DNAStringSet")
    {
        if (is.null(.y))
        {
            res = try(memes:runAme(
                .data,
                database=database
            ))
        }
        else
        {
            res = try(memes::runAme(
                .data, control=.y,
                database=database
            ))
        }
    }
    else {
       stop("Bad data type (not of class data.frame or DNAStringSet)")
    }

    if (is.null(res) | "try-error" %in% class(res))
    {
        return("No motifs found")
    }
    return(res)
}