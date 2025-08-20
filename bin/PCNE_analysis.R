#!/usr/bin/env Rscript

# Version: 1.0.0 compatible with PCNE v2.0.0

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  stop("Usage: Rscript PCNE_analysis.R <window_data.tsv> <chr_list.txt> <plasmid_list.txt> <enable_gc_flag> <loess_frac_str> <output.tsv> <generate_gc_plot_flag> <plot_flag> <agg_flag> <plasmid_fasta> <sample_name>", call. = FALSE)
}

window_data_file <- args[1]
chr_list_file <- args[2]
plasmid_list_file <- args[3]
enable_gc_flag <- as.logical(as.integer(args[4]))
loess_frac_str <- args[5]
output_file <- args[6]
generate_gc_plot_flag <- as.logical(as.integer(args[7]))
generate_plot_flag <- as.logical(as.integer(args[8]))
aggregate_flag <- as.logical(as.integer(args[9]))
plasmid_input_filename <- args[10]
sample_name <- args[11]

# --- Read Input Data ---
message("Reading windowed coverage data...")
window_data <- readr::read_tsv(window_data_file, show_col_types = FALSE) %>%
               filter(!is.na(depth), !is.na(gc), is.finite(depth), is.finite(gc), depth >= 0)
chr_names <- readLines(chr_list_file)
plasmid_names <- readLines(plasmid_list_file)

# Initialize the column to be used for final calculations
window_data$depth_to_use <- window_data$depth
final_norm_mode <- "Whole_chromosome"
gc_correction_applied <- FALSE
loess_frac <- NA
gc_model <- NULL

# --- Perform Unified LOESS GC Correction (if enabled) ---
if (enable_gc_flag) {
  message("GC Correction enabled. Building LOESS model...")

  baseline_windows_for_fit <- window_data %>% filter(contig %in% chr_names, depth > 0)

  if(nrow(baseline_windows_for_fit) < 50) {
      message("Warning: Too few baseline windows with coverage (<50). Skipping GC correction.")
  } else {
      # --- Automatic LOESS Fraction Selection with K-Fold Cross-Validation ---
      if (tolower(loess_frac_str) == "auto") {
          message("Selecting LOESS fraction via 5-fold cross-validation...")
          candidate_fracs <- seq(0.15, 0.75, by = 0.05) 
          k_folds <- 5
          
          # Create folds
          set.seed(123) # for reproducibility
          folds <- sample(cut(seq(1, nrow(baseline_windows_for_fit)), breaks = k_folds, labels = FALSE))
          
          avg_mse_per_frac <- sapply(candidate_fracs, function(span) {
              fold_mses <- sapply(1:k_folds, function(k) {
                  test_indices <- which(folds == k, arr.ind = TRUE)
                  train_data <- baseline_windows_for_fit[-test_indices, ]
                  test_data <- baseline_windows_for_fit[test_indices, ]
                  
                  model <- tryCatch(
                      loess(depth ~ gc, data = train_data, span = span, control = loess.control(surface = "direct")),
                      error = function(e) NULL
                  )
                  
                  if (is.null(model)) return(NA)
                  
                  predictions <- predict(model, newdata = test_data)
                  # Calculate Mean Squared Error for this fold
                  mean((test_data$depth - predictions)^2, na.rm = TRUE)
              })
              # Average the MSE across all k folds for this span
              mean(fold_mses, na.rm = TRUE)
          })
          
          best_frac_index <- which.min(avg_mse_per_frac)
          
          if (length(best_frac_index) == 0 || is.infinite(avg_mse_per_frac[best_frac_index])) {
              message("Warning: LOESS selection failed. Defaulting to 0.3")
              loess_frac <- 0.3
          } else {
              loess_frac <- candidate_fracs[best_frac_index]
              message(sprintf("Selected optimal LOESS fraction: %.2f (Avg. MSE=%.4f)", loess_frac, min(avg_mse_per_frac, na.rm=TRUE)))
          }
      } else {
          loess_frac <- as.numeric(loess_frac_str)
      }

      message(sprintf("Fitting LOESS model (span = %.2f) using %d baseline windows.", loess_frac, nrow(baseline_windows_for_fit)))
      
      gc_model <- loess(depth ~ gc, data = baseline_windows_for_fit, span = loess_frac,
                        control = loess.control(surface = "direct"))

      window_data$expected_depth <- predict(gc_model, newdata = window_data)
      window_data$expected_depth[is.na(window_data$expected_depth) | window_data$expected_depth <= 0] <- 0.01

      global_median_depth <- median(window_data$depth[window_data$depth > 0], na.rm = TRUE)
      
      message(sprintf("Global median depth: %.3f", global_median_depth))

      window_data$depth_to_use <- pmax(0, window_data$depth * (global_median_depth / window_data$expected_depth))
      
      final_norm_mode <- "Whole_chromosome_GC_Corrected"
      gc_correction_applied <- TRUE
      message("GC correction applied successfully.")
  }
}

# --- Aggregate Depths and Calculate PCN using MEDIAN ---
message("Aggregating depths...")

contig_summary <- window_data %>%
    group_by(contig) %>%
    summarise(
        length = sum(end - start, na.rm = TRUE),
        median_depth = median(depth_to_use, na.rm = TRUE),
        .groups = 'drop'
    )

baseline_windows_final <- window_data %>% filter(contig %in% chr_names)
baseline_depth <- median(baseline_windows_final$depth_to_use, na.rm = TRUE)

message(sprintf("Final baseline depth: %.3f", baseline_depth))

plasmid_raw_report <- contig_summary %>%
    filter(contig %in% plasmid_names) %>%
    mutate(
        baseline_median_depth = baseline_depth,
        normalization_mode = final_norm_mode,
        estimated_copy_number = ifelse(baseline_depth > 0, round(median_depth / baseline_depth, 2), NA_real_)
    ) %>%
    rename(plasmid_contig = contig)

# --- Handle Aggregation ---
if (aggregate_flag) {
    message("Aggregating results for all plasmid contigs...")
    agg_data <- plasmid_raw_report %>%
        summarise(
            total_length = sum(length, na.rm = TRUE),
            total_weighted_depth = sum(median_depth * length, na.rm = TRUE)
        )
    agg_median_depth <- if(agg_data$total_length > 0) agg_data$total_weighted_depth / agg_data$total_length else 0
    agg_plasmid_id <- if (nzchar(plasmid_input_filename)) sub("\\.[^.]*$", "", basename(plasmid_input_filename)) else "Aggregated_Plasmid"
    
    final_report <- tibble(
        sample = sample_name,
        plasmid_contig = agg_plasmid_id,
        length = as.integer(agg_data$total_length),
        median_depth = agg_median_depth,
        baseline_median_depth = baseline_depth,
        normalization_mode = final_norm_mode,
        estimated_copy_number = ifelse(baseline_depth > 0, round(agg_median_depth / baseline_depth,2), NA_real_)
    )
} else {
    final_report <- plasmid_raw_report %>%
      mutate(sample = sample_name) %>%
      select(sample, plasmid_contig, length, median_depth, baseline_median_depth, normalization_mode, estimated_copy_number) %>%
      arrange(desc(estimated_copy_number))
}

# --- Generate Plots ---
if (generate_plot_flag && nrow(final_report) > 0) {
    message("Generating final results plot...")
    output_prefix <- sub("_results\\.tsv$", "", output_file)
    plot_file <- paste0(output_prefix, "_plot.png")
    
    plot_data <- final_report %>%
        mutate(plot_label = factor(plasmid_contig, levels = rev(final_report$plasmid_contig)))
        
    p <- ggplot(plot_data, aes(x = plot_label, y = estimated_copy_number)) +
      geom_bar(stat = "identity", fill = "skyblue", color = "black", width=0.7) +
      coord_flip() +
      labs(
        title = "Plasmid Copy Numbers",
        subtitle = paste("Sample:", sample_name, "| Normalization:", final_norm_mode),
        x = "Plasmid",
        y = "Copy Number (Median-based)"
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, face="bold"),
        plot.subtitle = element_text(hjust = 0.5, size=9)
      ) +
      scale_y_continuous(expand = c(0, 0))
      
    ggsave(filename = plot_file, plot = p, width = 7, height = max(4, nrow(plot_data) * 0.25 + 1.5), limitsize = FALSE, dpi=300)
    message(paste("Plot saved successfully to:", plot_file))
}

# --- Generate GC diagnostic plot ---
if (gc_correction_applied && generate_gc_plot_flag && !is.null(gc_model)) {
    output_prefix <- sub("_results\\.tsv$", "", output_file)
    gc_plot_file <- paste0(output_prefix, "_gcplot.png")
    
    message(paste("Generating GC diagnostic plot to:", gc_plot_file))
    
    data_for_plot <- window_data %>% 
        filter(depth > 0) %>%
        mutate(source = if_else(contig %in% chr_names, "Baseline", "Plasmid"))
    data_for_plot$predicted_depth <- predict(gc_model, newdata = data_for_plot)
    
    p_gc <- ggplot(data_for_plot, aes(x = gc, y = depth)) +
      geom_point(aes(color = source), alpha = 0.3, size = 0.5) +
      scale_color_manual(values = c("Baseline" = "grey50", "Plasmid" = "dodgerblue")) +
      geom_line(aes(y = predicted_depth), color = "red", linewidth = 1) +
      scale_y_log10() +
      labs(
        title = "GC Content vs. Read Depth ",
        subtitle = paste("Sample:", sample_name, "| LOESS span =", round(loess_frac, 2)),
        x = "GC Fraction",
        y = "Mean Depth per Window (Log Scale)"
      ) +
      theme_classic(base_size = 11) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
    
    ggsave(filename = gc_plot_file, plot = p_gc, width = 7, height = 5, dpi = 300)
    message("GC diagnostic plot saved successfully.")
}

# --- Write Output ---
message(paste("Writing final report to:", output_file))
readr::write_tsv(final_report, output_file, na = "NA")

message("R script finished.")
