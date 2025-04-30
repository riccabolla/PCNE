#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  message("Usage: Rscript PCNE.R <plasmid_coverage.tsv> <denominator_depth> <norm_mode_string> <output.tsv> <generate_plot_flag> <aggregate_plasmid_flag> <plasmid_input_filename_or_empty>")
  message("  <plasmid_coverage.tsv>: Pre-filtered TSV file with columns: plasmid_contig, length, mean_depth")
  message("  <denominator_depth>: Single numeric value for baseline average depth")
  message("  <norm_mode_string>: String indicating normalization mode (e.g., 'Chromosome', 'BUSCO_SCG')")
  message("  <output.tsv>: Path for the final TSV report")
  message("  <generate_plot_flag>: 0 (no plot) or 1 (generate plot)")
  message("  <aggregate_plasmid_flag>: 0 (per-contig results) or 1 (aggregate results for all contigs in input)")
  message("  <plasmid_input_filename_or_empty>: Original plasmid input filename (-p value) if aggregating, otherwise \"\"")
  stop("Incorrect number of arguments supplied.", call. = FALSE)
}

plasmid_cov_file <- args[1]
denominator_depth_str <- args[2]
norm_mode_str <- args[3]
output_file <- args[4]
generate_plot_flag <- as.integer(args[5])
aggregate_flag <- as.logical(as.integer(args[6]))
plasmid_input_filename <- args[7]


denominator_depth <- as.numeric(denominator_depth_str)
if(is.na(denominator_depth)) {
    warning("Warning: Invalid baseline denominator depth value received. Setting to 0.", call.=FALSE)
    denominator_depth <- 0
}
message(paste("Using baseline mean depth:", round(denominator_depth, 2), "calculated via", norm_mode_str, "mode."))

plot_requested <- (generate_plot_flag == 1)


message(paste("Reading plasmid coverage data from:", plasmid_cov_file))
plasmid_data <- NULL
tryCatch({
    if (!file.exists(plasmid_cov_file)) { stop(paste("Plasmid coverage temp file not found:", plasmid_cov_file), call.=FALSE) }
    if (file.info(plasmid_cov_file)$size == 0) {
        warning(paste("Plasmid coverage temp file is empty:", plasmid_cov_file), call.=FALSE)
        plasmid_data <- tibble(plasmid_contig=character(), length=integer(), mean_depth=double())
    } else {
        plasmid_data <- readr::read_tsv(plasmid_cov_file, show_col_types = FALSE)
    }
    required_cols <- c("plasmid_contig", "length", "mean_depth")
    if(!all(required_cols %in% names(plasmid_data))) {
        print(paste("Available columns in plasmid data:", paste(names(plasmid_data), collapse=", ")))
        missing_cols <- setdiff(required_cols, names(plasmid_data))
        stop(paste("Plasmid coverage file missing required columns:", paste(missing_cols, collapse=", ")), call.=FALSE)
    }
    if(nrow(plasmid_data) == 0 && file.info(plasmid_cov_file)$size > 0) {
         warning("Warning: Plasmid coverage data file contains no valid data rows (only header?).", call. = FALSE)
    }
}, error = function(e){
    stop(paste("Error reading plasmid coverage file '", plasmid_cov_file, "': ", e$message), call.=FALSE)
})


message("Calculating copy numbers...")

output_columns <- c("plasmid_contig", "length", "mean_depth", "baseline_mean_depth",
                    "normalization_mode", "estimated_copy_number")

final_report <- data.frame(matrix(ncol = length(output_columns), nrow = 0,
                          dimnames = list(NULL, output_columns))) %>%
                mutate(
                    plasmid_contig=character(), length=integer(), mean_depth=double(),
                    baseline_mean_depth=double(), normalization_mode=character(),
                    estimated_copy_number=double()
                )

if (inherits(plasmid_data, "data.frame") && nrow(plasmid_data) > 0) {

    if (aggregate_flag) {
        message("Aggregating results for all plasmid contigs...")
        agg_data <- plasmid_data %>%
            summarise(
                total_length = sum(length, na.rm = TRUE),
                total_weighted_depth = sum(mean_depth * length, na.rm = TRUE),
                .groups = 'drop'
            )

        aggregated_mean_depth <- NA_real_
        estimated_copy_number <- NA_real_

        if (!is.na(agg_data$total_length) && agg_data$total_length > 0) {
            aggregated_mean_depth = agg_data$total_weighted_depth / agg_data$total_length
            estimated_copy_number = if_else(!is.na(denominator_depth) && denominator_depth > 0, aggregated_mean_depth / denominator_depth, NA_real_)

            agg_plasmid_id <- if (nzchar(plasmid_input_filename)) {
                                  sub("\\.[^.]*$", "", basename(plasmid_input_filename))
                              } else { "Aggregated_Plasmid" }

            final_report <- tibble(
                plasmid_contig = agg_plasmid_id,
                length = as.integer(agg_data$total_length),
                mean_depth = aggregated_mean_depth,
                baseline_mean_depth = denominator_depth,
                normalization_mode = norm_mode_str,
                estimated_copy_number = estimated_copy_number
            )
        } else { warning("Could not calculate aggregated metrics.", call. = FALSE) }

    } else {
        message("Calculating per-contig results...")
        final_report <- plasmid_data %>%
          dplyr::mutate(
            baseline_mean_depth = denominator_depth,
            normalization_mode = norm_mode_str,
            estimated_copy_number = if_else(!is.na(baseline_mean_depth) && baseline_mean_depth > 0, mean_depth / baseline_mean_depth, NA_real_)
          ) %>%
          dplyr::select( all_of(output_columns) ) %>% 
          dplyr::arrange(desc(estimated_copy_number))
    }
} else {
    warning("No valid plasmid coverage data found to calculate copy numbers.", call. = FALSE)
}


# --- Plot ---
if (plot_requested) {
    message("Generating plot...")
    suppressPackageStartupMessages(library(ggplot2))

    output_prefix <- sub("_results\\.tsv$", "", output_file)
    plot_file <- paste0(output_prefix, "_plot.png")

    if (nrow(final_report) > 0 && sum(!is.na(final_report$estimated_copy_number)) > 0) {
        plot_data <- final_report %>%
            mutate(
                plot_label = factor(plasmid_contig, levels = rev(final_report$plasmid_contig))
            )

        p <- ggplot(plot_data, aes(x = plot_label, y = estimated_copy_number)) +
          geom_bar(stat = "identity", fill = "skyblue", color = "black", width=0.7) +
          coord_flip() +
          labs(
            title = paste("Estimated Plasmid Copy Numbers"),
            x = if_else(aggregate_flag, "Aggregated Plasmid", "Plasmid Contig"),
            y = "Copy Number"
          ) +
          theme_minimal(base_size = 11) +
          theme(
              plot.title = element_text(hjust = 0.5, face="bold"),
              axis.text.y = element_text(size = if_else(aggregate_flag, 10, 8)),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
          )

        tryCatch({
            plot_height = max(4, nrow(plot_data) * 0.25 + 1.5)
            ggsave(filename = plot_file, plot = p, width = 7, height = plot_height, limitsize = FALSE, dpi=300)
            message(paste("Plot saved successfully to:", plot_file))
        }, error = function(e) {
            warning(paste("Failed to save plot:", e$message), call. = FALSE)
        })
    } else {
        warning("No data to plot.", call. = FALSE)
    }
}

# --- Write Output ---
message(paste("Writing final report to:", output_file))
final_report_rounded <- final_report %>%
    mutate(across(c(mean_depth, baseline_mean_depth, estimated_copy_number),
                  ~ round(., digits = 3)))

tryCatch({
    readr::write_tsv(final_report_rounded, output_file, na = "NA", col_names = TRUE)
    message("Output file written successfully.")
}, error = function(e) {
    stop(paste("Error writing output file '", output_file, "': ", e$message), call. = FALSE)
})

message("R script finished.")

