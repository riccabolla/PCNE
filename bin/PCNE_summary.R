#!/usr/bin/env Rscript

# Version: 1.0.0
# Summarizes results from multiple pcne runs.

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

# Read files
combined_data <- result_files %>%
  set_names() %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE), .id = "source_file") %>%
  mutate(sample = basename(sample),
    source_file = basename(source_file))

if (!("sample" %in% names(combined_data))) {
    stop("Error: The input TSV files do not contain the required 'sample' column. Please re-run pcne with a version that includes this feature.", call. = FALSE)
}

# --- Write Combined Table ---
output_summary_file <- "pcne_summary_all_results.tsv"
message(paste("Writing combined data to:", output_summary_file))
write_tsv(combined_data, output_summary_file)

# --- Generate Summary Plot ---
output_plot_file <- "pcne_summary_plot.png"
message(paste("Generating summary plot to:", output_plot_file))

# Create a summary plot
p <- ggplot(combined_data, aes(x = sample, y = estimated_copy_number, fill = plasmid_contig)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(
    title = "PCNE summary",
    x = "Sample",
    y = "Estimated Copy Number",
    fill = "Plasmid"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_y_continuous(expand = c(0, 0))

ggsave(filename = output_plot_file, plot = p, width = max(8, length(unique(combined_data$sample)) * 0.5), height = 6, dpi = 300, limitsize = FALSE)

message("Summary complete.")
