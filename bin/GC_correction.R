#!/usr/bin/env Rscript

# GC correction script for PCNE v1.0.0 
# Performs LOESS GC correction, calculates copy numbers, and plots.

# --- Load Libraries ---
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 10) {
  message("Usage: Rscript gc_correct_and_calculate.R <window_data.tsv> <baseline_list.txt> <plasmid_list.txt> <loess_frac> <output.tsv> <gc_plot_file_or_empty> <generate_plot_flag> <aggregate_plasmid_flag> <plasmid_input_filename_or_empty> <norm_mode_string>")
  message("  <window_data.tsv>: TSV with windowed contig, start, end, gc, depth")
  message("  <baseline_list.txt>: File containing list of contig names for baseline (Chr or BUSCO)")
  message("  <plasmid_list.txt>: File containing list of plasmid contig names")
  message("  <loess_frac>: LOESS span parameter (e.g., 0.3)")
  message("  <output.tsv>: Path for the final TSV report")
  message("  <gc_plot_file_or_empty>: Path to save GC diagnostic plot, or \"\"")
  message("  <generate_plot_flag>: 0 (no plot) or 1 (generate final plot)")
  message("  <aggregate_plasmid_flag>: 0 (per-contig) or 1 (aggregate)")
  message("  <plasmid_input_filename_or_empty>: Original plasmid input filename if aggregating")
  message("  <norm_mode_string>: String indicating normalization mode (e.g., Chromosome_GC_Corrected)")
  stop("Incorrect number of arguments supplied.", call. = FALSE)
}

window_data_file <- args[1]
baseline_list_file <- args[2]
plasmid_list_file <- args[3]
loess_frac_str <- args[4]
output_file <- args[5]
gc_plot_file <- args[6]
gc_plot_file <- trimws(gc_plot_file)
generate_plot_flag <- as.integer(args[7])
aggregate_flag <- as.logical(as.integer(args[8]))
plasmid_input_filename <- args[9]
norm_mode_str <- args[10]

# --- Inputs ---
loess_frac <- as.numeric(loess_frac_str)
if(is.na(loess_frac) || loess_frac <= 0 || loess_frac > 1) {
    warning(paste("Invalid LOESS fraction provided:", loess_frac_str, "- using default 0.3"), call.=FALSE)
    loess_frac <- 0.3
}
plot_requested <- (generate_plot_flag == 1)
gc_plot_requested <- nzchar(gc_plot_file) 

# --- Function to Read Name Lists ---
read_name_list <- function(filepath, list_type) {
  if (!file.exists(filepath) || file.access(filepath, 4) != 0) {
      stop(paste("Error:", list_type, "list file not found or not readable:", filepath), call. = FALSE)
  }
  tryCatch({
      names <- readLines(filepath)
      names <- trimws(names); names <- names[names != ""]
      if (length(names) == 0) { warning(paste("Warning: No names found in", list_type, "list file:", filepath), call. = FALSE) }
      return(unique(names))
  }, error = function(e) { stop(paste("Error reading", list_type, "list file '", filepath, "': ", e$message), call. = FALSE) })
}

# --- Read Lists ---
message(paste("Reading baseline contig names from:", baseline_list_file))
baseline_contig_names <- read_name_list(baseline_list_file, "Baseline")
message(paste("Found", length(baseline_contig_names), "unique baseline contig names."))

message(paste("Reading plasmid contig names from:", plasmid_list_file))
plasmid_names_list <- read_name_list(plasmid_list_file, "Plasmid")
message(paste("Found", length(plasmid_names_list), "unique plasmid contig names."))

# --- Read Window Data ---
message(paste("Reading windowed GC/depth data from:", window_data_file))
window_data <- NULL
tryCatch({
    if (!file.exists(window_data_file)) { stop(paste("Window data file not found:", window_data_file), call.=FALSE) }
    if (file.info(window_data_file)$size == 0) { stop(paste("Window data file is empty:", window_data_file), call.=FALSE) }

    window_data <- readr::read_tsv(window_data_file, show_col_types = FALSE)
    required_cols <- c("contig", "start", "end", "gc", "depth")
    if(!all(required_cols %in% names(window_data))) {
        stop(paste("Window data file missing required columns:", paste(setdiff(required_cols, names(window_data)), collapse=", ")), call.=FALSE)
    }
    window_data <- window_data %>% dplyr::mutate(length = end - start)
    window_data <- window_data %>% dplyr::filter(!is.na(depth), !is.na(gc), length > 0, is.finite(depth), is.finite(gc), depth >= 0)
    if(nrow(window_data) == 0) { stop("No valid data rows found in window data file after filtering NAs/zero length/neg depth.", call.=FALSE) }

}, error = function(e){ stop(paste("Error reading or processing window data file '", window_data_file, "': ", e$message), call.=FALSE) })
message("Windowed GC/depth data read successfully.")

# --- Perform GC Correction ---
message("Performing LOESS GC correction...")
baseline_windows <- window_data %>% dplyr::filter(contig %in% baseline_contig_names)

if(nrow(baseline_windows) < 10) {
    stop("Error: Too few valid baseline windows (<10) found for GC correction model fitting.", call.=FALSE)
}

message(paste("Fitting LOESS model (span =", loess_frac, ") using", nrow(baseline_windows), "baseline windows..."))
baseline_windows_for_fit <- baseline_windows %>% filter(depth > 0)
if(nrow(baseline_windows_for_fit) < 10) {
     stop("Error: Too few baseline windows with non-zero depth (<10) found for GC correction model fitting.", call.=FALSE)
}

gc_model <- tryCatch({
    stats::loess(depth ~ gc, data = baseline_windows_for_fit, span = loess_frac,
                 control = loess.control(surface = "direct"))
}, error = function(e) {
    stop(paste("Error fitting LOESS model:", e$message), call.=FALSE)
})

message("Predicting expected depth for all windows...")
window_data <- window_data %>%
    dplyr::mutate(
      expected_depth = suppressWarnings(predict(gc_model, newdata = gc)),
      expected_depth = if_else(is.na(expected_depth) | expected_depth <= 0, 0.01, expected_depth)
    )

global_median_baseline_depth <- median(baseline_windows_for_fit$depth, na.rm = TRUE)
if(is.na(global_median_baseline_depth) || global_median_baseline_depth <= 0) {
    warning("Warning: Could not calculate valid global median baseline depth. Correction might be unreliable.", call.=FALSE)
    global_median_baseline_depth <- 1
}
message(paste("Global median baseline depth (from fit data):", round(global_median_baseline_depth, 2)))

message("Calculating corrected depth for all windows...")
window_data <- window_data %>%
    dplyr::mutate(
      corrected_depth = pmax(0, depth * (global_median_baseline_depth / expected_depth))
    )

# --- Aggregate Corrected Depths ---
message("Aggregating corrected depths per contig...")
contig_corrected_summary <- window_data %>%
    dplyr::group_by(contig) %>%
    dplyr::summarise(
        total_length = sum(length, na.rm = TRUE),
        total_weighted_corrected_depth = sum(corrected_depth * length, na.rm = TRUE),
        .groups = 'drop'
    ) %>%
    dplyr::mutate(
        corrected_mean_depth = if_else(total_length > 0, total_weighted_corrected_depth / total_length, 0)
    )

baseline_summary <- contig_corrected_summary %>%
    dplyr::filter(contig %in% baseline_contig_names)

denominator_depth_corrected <- 0
if(nrow(baseline_summary) > 0 && sum(baseline_summary$total_length, na.rm=TRUE) > 0) {
   denominator_depth_corrected <- sum(baseline_summary$total_weighted_corrected_depth, na.rm=TRUE) / sum(baseline_summary$total_length, na.rm=TRUE)
} else {
   warning("Warning: Could not calculate corrected baseline denominator depth (no baseline contigs found or zero total length).", call. = FALSE)
}
message(paste("Final baseline mean depth (GC Corrected):", round(denominator_depth_corrected, 2)))

plasmid_corrected_data <- contig_corrected_summary %>%
    dplyr::filter(contig %in% plasmid_names_list) %>%
    dplyr::select(plasmid_contig = contig, length = total_length, mean_depth = corrected_mean_depth)

# --- Calculate Copy Number  ---
message("Calculating copy numbers using corrected depths (Poisson CIs removed)...")

# Define output columns 
output_columns <- c("plasmid_contig", "length", "mean_depth", "baseline_mean_depth",
                    "normalization_mode", "estimated_copy_number")

# Initialize report
final_report <- data.frame(
    plasmid_contig = character(),
    length = integer(),
    mean_depth = double(),
    baseline_mean_depth = double(),
    normalization_mode = character(),
    estimated_copy_number = double(),
    stringsAsFactors = FALSE
)


if (inherits(plasmid_corrected_data, "data.frame") && nrow(plasmid_corrected_data) > 0) {

    if (aggregate_flag) {
        message("Aggregating GC-corrected results for all plasmid contigs...")
        agg_data <- plasmid_corrected_data %>%
            summarise(
                total_length = sum(length, na.rm = TRUE),
                total_weighted_depth = sum(mean_depth * length, na.rm = TRUE),
                .groups = 'drop'
            )
        aggregated_mean_depth <- NA_real_
        estimated_copy_number <- NA_real_
        if (!is.na(agg_data$total_length) && agg_data$total_length > 0) {
            aggregated_mean_depth = agg_data$total_weighted_depth / agg_data$total_length
            estimated_copy_number = if_else(!is.na(denominator_depth_corrected) & denominator_depth_corrected > 0, aggregated_mean_depth / denominator_depth_corrected, NA_real_)
            agg_plasmid_id <- if (nzchar(plasmid_input_filename)) { sub("\\.[^.]*$", "", basename(plasmid_input_filename)) } else { "Aggregated_Plasmid" }            
            final_report <- tibble(
                plasmid_contig = agg_plasmid_id,
                length = as.integer(agg_data$total_length),
                mean_depth = aggregated_mean_depth,
                baseline_mean_depth = denominator_depth_corrected,
                normalization_mode = norm_mode_str,
                estimated_copy_number = estimated_copy_number
            )
        } else { warning("Could not calculate aggregated metrics.", call. = FALSE) }

    } else {
        message("Calculating per-contig results using corrected depths...")
        final_report <- plasmid_corrected_data %>%
          dplyr::mutate( 
            baseline_mean_depth = denominator_depth_corrected, 
            normalization_mode = norm_mode_str, 
            estimated_copy_number = if_else(!is.na(baseline_mean_depth) & baseline_mean_depth > 0, mean_depth / baseline_mean_depth, NA_real_)
          ) %>%
          dplyr::select( all_of(output_columns) ) %>%
          dplyr::arrange(desc(estimated_copy_number))
    }
} else {
    warning("No valid plasmid coverage data found.", call. = FALSE)
}


# --- Generate Plots ---
if (plot_requested || gc_plot_requested) {
    message("Loading ggplot2 for plotting...")
    suppressPackageStartupMessages(library(ggplot2))
}

if (gc_plot_requested) {
    message(paste("Attempting to generate GC diagnostic plot to file:", gc_plot_file))
    message(paste("R's current working directory is:", getwd())) 
    tryCatch({
        plot_gc_data <- window_data %>%
            mutate(Baseline = if_else(contig %in% baseline_contig_names, "Baseline", "Other"))

        p_gc <- ggplot(plot_gc_data, aes(x = gc, y = depth)) +
          geom_point(aes(color = Baseline), alpha = 0.3, size=0.5) +
          geom_line(aes(y = expected_depth), color = "red", linewidth = 0.8) +
          scale_color_manual(values=c("Legend"="blue", "Other"="purple")) +
          scale_y_log10(limits=c(0.1, NA), oob = scales::squish) +
          labs(
            title = "GC Content",
            subtitle = paste("LOESS=", loess_frac),
            x = "GC Fraction",
            y = "Mean Depth per Window (Log Scale)"
          ) +
          theme_classic(base_size = 11) +
          theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size=9))
        if (nzchar(gc_plot_file)) {
            ggsave(filename = gc_plot_file, plot = p_gc, width = 7, height = 5, dpi = 300)
            message(paste("GC diagnostic plot successfully saved to:", gc_plot_file))
        } else {
            warning("GC plot requested, but the provided filename was empty. Plot not saved.", call. = FALSE)
        }
    }, error = function(e) {
        warning(paste("Failed to generate GC diagnostic plot:", e$message, "(filename was '", gc_plot_file, "')"), call. = FALSE)
    })
}

if (plot_requested) {
    message("Generating final results plot...")
    output_prefix <- sub("_results\\.tsv$", "", output_file)
    plot_file <- paste0(output_prefix, "_plot.png") 
    if (inherits(final_report, "data.frame") && nrow(final_report) > 0 && sum(!is.na(final_report$estimated_copy_number)) > 0) {
        plot_data <- final_report %>%
            mutate(
                plot_label = factor(plasmid_contig, levels = rev(unique(final_report$plasmid_contig[order(final_report$estimated_copy_number)])))
            )
        p <- ggplot(plot_data, aes(x = plot_label, y = estimated_copy_number)) +
          geom_bar(stat = "identity", fill = "skyblue", color = "black", width=0.7) +
          coord_flip() +
          labs(
            title = paste("Plasmid Copy Numbers"),
            subtitle = paste("Normalization Mode:", norm_mode_str),
            x = "Plasmid",
            y = "Copy Number"
          ) +
          scale_y_continuous(expand = c(0, 0)) +
          theme_classic(base_size = 11) +
          theme(
            plot.title = element_text(hjust = 0.5, face="bold"),
            plot.subtitle = element_text(hjust = 0.5, size=9),
            axis.text.y = element_text(size = if_else(nrow(plot_data) > 20, 6, if_else(nrow(plot_data) > 10, 8, 10))),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
          )

        tryCatch({
            plot_height = max(4, nrow(plot_data) * 0.25 + 1.5)
            if (nzchar(plot_file)) { 
                ggsave(filename = plot_file, plot = p, width = 7, height = plot_height, limitsize = FALSE, dpi=300)
                message(paste("Plot saved successfully to:", plot_file))
            } else {
                 warning("Final results plot filename was empty. Plot not saved.", call. = FALSE)
            }
        }, error = function(e) {
            warning(paste("Failed to save plot:", e$message), call. = FALSE)
        })
    } else {
        warning("No data to plot for final results.", call. = FALSE)
    }
}

# --- Write Output ---
message(paste("Writing final report to:", output_file))

if (inherits(final_report, "data.frame") && nrow(final_report) > 0) {
    numeric_cols_to_round <- intersect(names(final_report), c("mean_depth", "baseline_mean_depth", "estimated_copy_number"))
    
    final_report_rounded <- final_report %>%
        mutate(across(all_of(numeric_cols_to_round), ~ round(., digits = 3)))

    tryCatch({
        readr::write_tsv(final_report_rounded, output_file, na = "NA", col_names = TRUE)
        message("Output file written successfully.")
    }, error = function(e) {
        stop(paste("Error writing output file '", output_file, "': ", e$message), call. = FALSE)
    })
} else {
     message("Final report is empty.")
     empty_df_for_output <- tibble(
        plasmid_contig = character(),
        length = integer(),
        mean_depth = double(),
        baseline_mean_depth = double(),
        normalization_mode = character(),
        estimated_copy_number = double()
     )
     tryCatch({
        readr::write_tsv(empty_df_for_output, output_file, na = "NA", col_names = TRUE)
        message("Empty output written successfully.")
    }, error = function(e) {
        stop(paste("Error writing empty output file '", output_file, "': ", e$message), call. = FALSE)
    })
}

message("R script finished.")

