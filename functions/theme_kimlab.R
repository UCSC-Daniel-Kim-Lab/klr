txt.mm_to_pts <- function(ratio){ base_size * ratio * (25.4 / 72.27) }
base_size = 8

te_aware_annotation_levels <- c('gencode_coding',
                                 'gencode_lncRNA',
                                 'LINE',
                                 'SINE',
                                 'LTR',
                                 'DNA',
                                 'Simple_repeat',
                                 'other')
te_color_pal <- c('grey', RColorBrewer::brewer.pal(n = 7, 'Dark2'))
names(te_color_pal) = te_aware_annotation_levels

te_alone_annotation_levels <- c('LINE',
                                'SINE',
                                'LTR',
                                'DNA',
                                'Simple_repeat',
                                'rRNA',
                                'tRNA')
te_alone_color_pal <- RColorBrewer::brewer.pal(n = 8, 'Dark2')[2:8]
names(te_alone_color_pal) <- te_alone_annotation_levels

#' The Official GGplot theme of the Daniel H. Kim Lab, UC Santa Cruz.
#' 
#' @export
theme_kimlab <- function()
{
    ggplot2::theme_set(
        ggthemes::theme_foundation(
            base_size = base_size,
            base_family = 'Helvetica'
        )+
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                size = base_size,
                face = 'bold'
            ),
            panel.background = ggplot2::element_rect(colour = NA),
            plot.background = ggplot2::element_rect(colour = NA),
            panel.border = ggplot2::element_rect(colour = NA), 
            axis.line = ggplot2::element_line(), 
            axis.line.x = NULL, 
            axis.line.y = NULL, 
            axis.text = ggplot2::element_text(size = rel(0.95)), 
            axis.text.x = ggplot2::element_text(
                margin = ggplot2::margin(t = 0.8 * base_size/4)), 
                axis.text.x.top = ggplot2::element_text(
                    margin = ggplot2::margin(b = 0.8 * base_size/4), 
                    vjust = 0
                ), 
                axis.text.y = ggplot2::element_text(
                    margin = ggplot2::margin(r = 0.5 * base_size/4), 
                    hjust = 1
                ),
                axis.text.y.right = ggplot2::element_text(
                    margin = ggplot2::margin(l = 0.5 * base_size/4), 
                    hjust = 0
                ), 
                axis.ticks = ggplot2::element_line(), 
                axis.ticks.length = ggplot2::unit(base_size/2.5, "pt"), axis.ticks.length.x = NULL, 
                axis.ticks.length.x.top = NULL, axis.ticks.length.x.bottom = NULL, 
                axis.ticks.length.y = NULL, axis.ticks.length.y.left = NULL, 
                axis.ticks.length.y.right = NULL,
                strip.text = ggplot2::element_text(size = rel(0.8), face = 'bold'),
                strip.background = ggplot2::element_blank(),
                legend.key.size= ggplot2::unit(0.08, "in"),
                legend.spacing = ggplot2::unit(0, "in"),
                legend.key = ggplot2::element_rect(colour = NA),
                legend.title = ggplot2::element_text(face="italic"),
                legend.text = ggplot2::element_text(face = 'bold'),
                legend.justification = c("right", "top"),
                legend.box.just = "right",
                legend.margin = ggplot2::margin(6, 6, 6, 6),
                plot.tag = ggplot2::element_text(
                    size = base_size,
                    face = 'bold'
                ),
                plot.margin = ggplot2::margin(0.04, 0.04, 0.04, 0.04, unit = "in"),
                legend.box.spacing = ggplot2::unit(-0.02, 'in'),
                panel.grid.major = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank()
        )
    )
}