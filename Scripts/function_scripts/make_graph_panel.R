library(cowplot)
library(magick)
library(ggplot2)
library(grid)

make_graph_panel <- function(png_files) {
  # Read the PNG files
  plots <- lapply(png_files, image_read)
  # Convert to grobs
  plot_grobs <- lapply(plots, function(x) {
    rasterGrob(as.raster(x))
  })
  
  # Create the main 2x2 grid
  left_panel <- plot_grid(
    plot_grobs[[1]], plot_grobs[[2]],
    plot_grobs[[3]], plot_grobs[[4]],
    labels = c("B", "C", "D", "E"),
    label_size = 20,
    ncol = 2
  )
  
  # Create the right panel with plot E
  right_panel <- plot_grid(
    NULL,
    plot_grobs[[5]],
    NULL,
    labels = c("","A", ""),
    label_size = 20,
    ncol = 1,
    rel_heights = c(0.2, 0.6, 0.2)
  )
  
  # Combine the panels with more space between them
  panel <- plot_grid(
    right_panel, left_panel,
    ncol = 2,
    rel_widths = c(0.7, 1.8),
    align = 'hv',
    axis = 'l',
    scale = c(1, 0.95)
    # spacing = 0.00
  )
  
  # Add border to the final panel with padding
  panel <- ggdraw(panel) + 
    theme(
      plot.background = element_rect(color = "black", fill = "white", linewidth = 1),
      plot.margin = margin(10, 10, 10, 10, "points")
    )


  return(panel)
}

# Example usage:
png_files <- c(
"dataset_specific_outputs/cecal/figures/cut_off_selection/phenotype/cecal_phenotype_rules_v_cutoffs_cosine.png",
"dataset_specific_outputs/cecal/figures/cut_off_selection/phenotype/cecal_phenotype_rules_v_cutoffs_jaccard.png",
"dataset_specific_outputs/cecal/figures/cut_off_selection/phenotype/cecal_phenotype_rules_v_cutoffs_kulczynski.png", 
"dataset_specific_outputs/cecal/figures/cut_off_selection/phenotype/cecal_phenotype_rules_v_cutoffs_support.png",
"dataset_specific_outputs/cecal/figures/rarefaction_curves/phenotype/cecal_phenotype_unique_classes_rules.png")


panel <- make_graph_panel(png_files)
save_plot("combined_outputs/combined_figures_outputs/Fig1_in_manuscript.png", panel, base_height = 10, base_width = 15)


print("DONE")
