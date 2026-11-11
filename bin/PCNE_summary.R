#!/usr/bin/env Rscript

# Version: 1.2.0 - Adaptive Faceted & Top-N Summary Plots

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(purrr))

message("Searching for pcne result files (*_results.tsv) in the current directory...")
result_files <- list.files(path = ".", pattern = "*_results\\.tsv$", full.names = TRUE)

if (length(result_files) == 0) {
  stop("No pcne result files (*_results.tsv) found in the current directory. Nothing to summarize.", call. = FALSE)
}

message(paste("Found", length(result_files), "result files. Combining..."))

combined_data <- result_files %>%
  set_names() %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE), .id = "source_file") %>%
  mutate(
    sample = basename(sample),
    source_file = basename(source_file)
  )

if (!("sample" %in% names(combined_data))) {
  stop("Error: The input TSV files do not contain the required 'sample' column.", call. = FALSE)
}

# --- Write Combined Table ---
output_summary_file <- "pcne_summary_all_results.tsv"
message(paste("Writing combined data to:", output_summary_file))
write_tsv(combined_data, output_summary_file)

# --- Count dimensions ---
n_samples  <- n_distinct(combined_data$sample)
n_plasmids <- n_distinct(combined_data$plasmid_contig)
n_rows     <- nrow(combined_data)

message(sprintf("Dataset contains %d samples, %d unique plasmid names (%d total entries).", 
                n_samples, n_plasmids, n_rows))

output_plot_file <- "pcne_summary_plot.png"

if (n_rows <= 20) {
  # --- SMALL DATASET: Standard Dodged Bar Plot ---
  message("<= 20 entries detected. Generating standard bar plot...")
  
  p <- ggplot(combined_data, aes(x = sample, y = estimated_copy_number, fill = plasmid_contig)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(
      title = "PCNE Summary",
      x = "Sample",
      y = "Estimated Copy Number",
      fill = "Plasmid"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  ggsave(output_plot_file, plot = p, width = max(8, n_samples * 0.6), height = 6, dpi = 300, limitsize = FALSE)
  
} else {
  # --- LARGE DATASET: Faceted Grid by Sample (No empty white space!) ---
  message("> 20 entries detected. Generating multi-panel faceted plot by sample...")
  
  # Order plasmids by copy number within each sample
  plot_data <- combined_data %>%
    mutate(unique_id = paste0(sample, "___", plasmid_contig)) %>%
    arrange(sample, estimated_copy_number) %>%
    mutate(unique_id = factor(unique_id, levels = unique(unique_id)))
  
  # Calculate grid dimensions (e.g., 4 or 5 columns wide)
  n_cols <- min(5, max(2, ceiling(sqrt(n_samples))))
  n_grid_rows <- ceiling(n_samples / n_cols)
  
  p_facet <- ggplot(plot_data, aes(x = unique_id, y = estimated_copy_number, fill = log10(estimated_copy_number + 1))) +
    geom_col(color = "black", width = 0.75) +
    geom_text(aes(label = round(estimated_copy_number, 1)), hjust = -0.15, size = 2.8) +
    coord_flip() +
    facet_wrap(~ sample, scales = "free_y", ncol = n_cols) +
    scale_x_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    scale_fill_viridis_c(option = "plasma", guide = "none") +
    labs(
      title = "PCNE Summary by Sample",
      subtitle = sprintf("Showing all %d plasmids across %d samples", n_rows, n_samples),
      x = "Plasmid Contig",
      y = "Estimated Copy Number"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid.minor = element_blank()
    )
  
  plot_width  <- n_cols * 4.2
  plot_height <- max(5, n_grid_rows * 2.2)
  
  ggsave(output_plot_file, plot = p_facet, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
  message(paste("Saved full faceted plot to:", output_plot_file))
  
  # --- Save a compact Top-30 Highest Copy Number Plot ---
  top_n_val <- 30
  if (n_rows > top_n_val) {
    top_plot_file <- "pcne_summary_top30_plot.png"
    message(paste("Also generating compact Top-30 plot to:", top_plot_file))
    
    top_data <- combined_data %>%
      arrange(desc(estimated_copy_number)) %>%
      slice_head(n = top_n_val) %>%
      mutate(label = paste0(sample, " | ", plasmid_contig)) %>%
      mutate(label = factor(label, levels = rev(unique(label))))
    
    p_top <- ggplot(top_data, aes(x = label, y = estimated_copy_number, fill = estimated_copy_number)) +
      geom_col(color = "black", width = 0.75) +
      geom_text(aes(label = round(estimated_copy_number, 1)), hjust = -0.15, size = 3) +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      scale_fill_viridis_c(name = "Copy Number", option = "plasma") +
      labs(
        title = sprintf("Top %d Highest Copy Number Plasmids", top_n_val),
        subtitle = sprintf("Filtered from %d total plasmids across %d samples", n_rows, n_samples),
        x = "Sample | Plasmid",
        y = "Estimated Copy Number"
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9)
      )
    
    ggsave(top_plot_file, plot = p_top, width = 9, height = 8, dpi = 300)
  }
}

message("Summary complete.")
