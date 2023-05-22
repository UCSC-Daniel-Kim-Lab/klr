#' Save a ggplot figure as a PDF
#' 
#' @param filename  The name of the file in which to save the panel.
#' @param figure_dir [Optional] The name of the directory in which to save the figure. Will be
#'  appended to the front of filename.
#' @param plot The GGPlot object to save.
#' @param annotation [Default=True] Whether or not to add the figure dimensions to the filename.
#' @param width The width of the panel to save.
#' @param height The height of the panel to save.
#' @param units The units to determine panel size.
#' @param device The device used to store the saved file.
#' @export
save_panel <- function(
    filename, figure_dir, plot=ggplot2::last_plot(), annotation=TRUE,
    width=90, height=180, units='mm', device=grDevices::cairo_pdf
)
{
    out_name = basename(filename)
    out_dir = substr(filename, 0, nchar(filename)-nchar(out_name))
    if (!is.null(figure_dir))
    {
        out_dir = paste0(figure_dir, "/", filename)
    }
    
    if (!dir.exists(out_dir)) dir.create(out_dir)

    ggsave(
        filename, plot,
        width=width, height=height, units=units,
        device=device
    )
}