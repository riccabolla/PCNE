#!/usr/bin/env Rscript

# R script for PCNE v0.1.0 / v0.1.1
# --- Load Libraries ---
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  message("Usage: Rscript PCNE.R <coverage_summary.txt> <chr_name_list.txt> <pls_name_list.txt> <output.tsv>")
  message("  <coverage_summary.txt>: Output file from 'samtools coverage'")
  message("  <chr_name_list.txt>: File containing list of chromosome contig names (one per line)")
  message("  <pls_name_list.txt>: File containing list of plasmid contig names (one per line)")
  message("  <output.tsv>: Path for the final TSV report")
  stop("Incorrect number of arguments supplied.", call. = FALSE)
}
coverage_file <- args[1]
chr_list_file <- args[2]
pls_list_file <- args[3]
output_file <- args[4]

# --- Function to Read Name Lists ---
read_name_list <- function(filepath, list_type) {
  if (!file.exists(filepath) || file.access(filepath, 4) != 0) { 
      stop(paste("Error:", list_type, "list file not found or not readable:", filepath), call. = FALSE)
  }
  tryCatch({
      names <- readLines(filepath)
      names <- trimws(names)
      names <- names[names != ""]
      if (length(names) == 0) {
          warning(paste("Warning: No names found in", list_type, "list file:", filepath), call. = FALSE)
      }
      return(unique(names))
  }, error = function(e) {
      stop(paste("Error reading", list_type, "list file '", filepath, "': ", e$message), call. = FALSE)
  })
}

# --- Read Name Lists ---
message(paste("Reading chromosome contig names from:", chr_list_file))
chr_names_list <- read_name_list(chr_list_file, "Chromosome")
message(paste("Found", length(chr_names_list), "unique chromosome contig names."))

message(paste("Reading plasmid contig names from:", pls_list_file))
plasmid_names_list <- read_name_list(pls_list_file, "Plasmid")
message(paste("Found", length(plasmid_names_list), "unique plasmid contig names."))

# --- Check for Overlap ---
overlap <- intersect(chr_names_list, plasmid_names_list)
if (length(overlap) > 0) {
    warning(paste("Warning: The following contigs are listed as both chromosome and plasmid:", paste(overlap, collapse=", ")), call. = FALSE)
    # Prioritize plasmid status by removing overlapping contigs from the chromosome list for normalization
    chr_names_list <- setdiff(chr_names_list, overlap)
    message(paste("Overlapping contigs will be treated as plasmids. Updated chromosome list size for normalization:", length(chr_names_list)))
}

# --- Read and Process Coverage Data ---
message(paste("Reading coverage file:", coverage_file))
coverage_data <- NULL # Initialize
tryCatch({
  # Define the exact column names based on the expected header from samtools coverage
  expected_colnames <- c("#rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")
  coverage_data <- readr::read_tsv(
      coverage_file,
      skip = 1,                   
      col_names = expected_colnames, 
      show_col_types = FALSE
  ) %>%
    {required_cols <- c("#rname", "startpos", "endpos", "meandepth")
      if (!all(required_cols %in% names(.))) {
        print(paste("Assigned column names:", paste(names(.), collapse=", ")))
        missing_cols <- setdiff(required_cols, names(.))
        stop(paste("Coverage file missing required columns after manual assignment:", paste(missing_cols, collapse=", ")), call. = FALSE)
      }
      . 
    } %>%
    dplyr::rename(contig_name = `#rname`) %>%
    dplyr::mutate(length = endpos - startpos,
                  across(c(startpos, endpos, length, meandepth), as.numeric)) %>%
    dplyr::filter(length > 0) %>%
    { if(nrow(.) == 0) stop("No valid contigs found in coverage data after filtering.", call. = FALSE); . }

}, error = function(e) {
  stop(paste("Error reading or processing coverage file '", coverage_file, "': ", e$message), call. = FALSE)
})
message("Coverage data read and processed successfully.")

# --- Calculate Weighted Chromosome Coverage (Denominator) ---
message("Calculating chromosome coverage...")
chromosome_data <- coverage_data %>%
  dplyr::filter(contig_name %in% chr_names_list)

denominator_avg_depth <- 0 
chromosome_summary <- list(contig_count = 0, total_length = 0) 

if (nrow(chromosome_data) == 0) {
    warning("Warning: No chromosome contigs from the list found in coverage data. Chromosome average depth will be 0.", call. = FALSE)
} else {
    chromosome_summary <- chromosome_data %>%
      dplyr::summarise(
        total_length = sum(length, na.rm = TRUE),
        total_weighted_depth = sum(meandepth * length, na.rm = TRUE),
        contig_count = n(),
        .groups = 'drop'
      )

    message(sprintf("Chromosome Summary (for Normalization): %d contigs processed, Total Length %d, Weighted Avg Depth Calculation...",
                chromosome_summary$contig_count %||% 0,
                chromosome_summary$total_length %||% 0
                ))

    if (!is.na(chromosome_summary$total_length) && chromosome_summary$total_length > 0) {
        denominator_avg_depth <- chromosome_summary$total_weighted_depth / chromosome_summary$total_length
        message(sprintf(" -> Weighted Avg Depth = %.2f", denominator_avg_depth))
    } else {
      warning("Total chromosome length calculated is zero, negative, or NA. Chromosome average depth will be 0.", call. = FALSE)
      denominator_avg_depth <- 0 
    }
}

baseline_mean_depth <- denominator_avg_depth
normalization_mode <- "Chromosome" 


# --- Calculate Copy Number for Each Plasmid ---
message("Calculating plasmid copy numbers...")
all_plasmid_coverage_data <- coverage_data %>%
  dplyr::filter(contig_name %in% plasmid_names_list)

results_list <- list()

if (length(plasmid_names_list) == 0) {
    warning("Plasmid name list is empty. No copy numbers to calculate.", call. = FALSE)
} else {
    if (nrow(all_plasmid_coverage_data) == 0) {
         warning("Warning: None of the plasmid names found in the list were present in the coverage summary file.", call. = FALSE)
    }
    for (p_name in plasmid_names_list) {
      plasmid_row <- all_plasmid_coverage_data %>% dplyr::filter(contig_name == p_name)

      if (nrow(plasmid_row) == 1) {
        p_depth <- plasmid_row$meandepth[1]
        p_length <- plasmid_row$length[1]
        copy_number <- if (!is.na(baseline_mean_depth) && baseline_mean_depth > 0) {
          p_depth / baseline_mean_depth
        } else {
          NA_real_ 
        }
        results_list[[p_name]] <- data.frame(
          plasmid_contig = p_name,
          length = p_length,
          mean_depth = p_depth,
          baseline_mean_depth = baseline_mean_depth, 
          normalization_mode = normalization_mode, 
          estimated_copy_number = copy_number
        )
      } else {
         warning(paste("Plasmid contig '", p_name, "' from list not found in coverage summary or had multiple entries. Skipping."), call. = FALSE)
          results_list[[p_name]] <- data.frame(
            plasmid_contig = p_name, length = NA_integer_, mean_depth = NA_real_,
            baseline_mean_depth = baseline_mean_depth,
            normalization_mode = normalization_mode,
            estimated_copy_number = NA_real_
          )
      }
    }
}

# --- Combine results and Write Output ---
output_columns <- c("plasmid_contig", "length", "mean_depth", "baseline_mean_depth",
                    "normalization_mode", "estimated_copy_number")

if (length(results_list) > 0) {
  final_report <- dplyr::bind_rows(results_list) %>%
      dplyr::select(all_of(output_columns)) %>%
      dplyr::arrange(desc(estimated_copy_number))

} else {
  warning("No plasmid copy numbers were calculated. Output file will be empty.", call. = FALSE)
  final_report <- data.frame(
      plasmid_contig=character(),
      length=integer(),
      mean_depth=double(),
      baseline_mean_depth=double(),
      normalization_mode=character(),
      estimated_copy_number=double()
     )
  names(final_report) <- output_columns
}

message(paste("Writing final report to:", output_file))
tryCatch({
    readr::write_tsv(final_report, output_file, na = "NA", col_names = TRUE)
    message("Output file written successfully.")
}, error = function(e) {
    stop(paste("Error writing output file '", output_file, "': ", e$message), call. = FALSE)
})

message("R script finished.")