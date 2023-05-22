#' Isolate the sequences in a set that contain motifs from a set list. 
#' 
#' @param .data A dataframe with a column of motif lists, or a motif list itself.
#' @param .obs The list of sequences to scan, used if a vector of motifs is provided.
#' @param motif_col The column of the data frame holding the motif lists.
#' @param seqs_col The column of the data frame holding the DNAStringSets of sequences to scan.
#' @export
source_motifs <- function(.data, .obs, motif_col="motifs", seqs_col="obs")
    if ("data.frame" %in% class(.data))
    {
        res = apply(.data, MARGIN=1, FUN=function(row)
        {
            temp_db = universalmotifs::write_meme(row[[motif_col]], "temp.meme", overwrite=F)
            temp = try(memes::runFimo(
                row[[seqs_col]], "temp.meme"
            ))
            file.remove("temp.meme")
            if (!("data.frame" %in% class(temp)))
            {
                temp = "No sequences matched"
            }
            return(temp)
        })
        return(res)
    } else if (!missing(.y))
    {
        temp_db = universalmotifs::write_meme(.data, "temp.meme", overwrite=F)
        temp = try(memes::runFimo(
            .obs, "temp.meme"
        ))
        file.remove("temp.meme")
        if (!("data.frame" %in% class(temp)))
        {
            temp = "No sequences matched"
        }
        return(temp)
    } else
    {
       stop("Invalid input data")
    }
}