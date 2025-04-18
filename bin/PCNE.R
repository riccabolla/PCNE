#!/usr/bin/env Rscript

# Script to calculate plasmid copy number
# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  message("Usage: Rscript calculate_multi_plasmid_cn.R <coverage_summary.txt> <chr_name_list.txt> <pls_name_list.txt> <output.tsv>")
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
  if (!file.exists(filepath)) {
      stop(paste("Error:", list_type, "list file not found:", filepath), call. = FALSE)
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
    chr_names_list <- setdiff(chr_names_list, overlap)
    message(paste("Overlapping contigs will be treated as plasmids. Updated chromosome list size:", length(chr_names_list)))
}

# --- Read and Process Coverage Data ---
message(paste("Reading coverage file:", coverage_file))
tryCatch({
  # Load libraries needed for this block
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(dplyr))

  # Define the exact column names based on the expected header
  expected_colnames <- c("#rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq")

  # Read the file using skip=1 and providing col_names manually
  coverage_data <- readr::read_tsv(
      coverage_file,
      skip = 1,                   
      col_names = expected_colnames,
      show_col_types = FALSE
  ) %>%
    {
      # Check if the required columns were assigned correctly
      required_cols <- c("#rname", "startpos", "endpos", "meandepth")

      if (!all(required_cols %in% names(.))) {
        print(paste("Assigned column names:", paste(names(.), collapse=", ")))
        missing_cols <- setdiff(required_cols, names(.))
        stop(paste("Coverage file missing required columns after manual assignment:", paste(missing_cols, collapse=", ")), call. = FALSE)
      }
      . # Pass data frame through if check passes
    } %>%
    # Rename using backticks
    dplyr::rename(contig_name = `#rname`) %>%
    # Mutate and filter using standard names
    dplyr::mutate(length = endpos - startpos,
                  across(c(startpos, endpos, length, meandepth), as.numeric)) %>%
    dplyr::filter(length > 0) %>%
    { if(nrow(.) == 0) stop("No valid contigs found in coverage data after filtering.", call. = FALSE); . }

}, error = function(e) {
  stop(paste("Error reading or processing coverage file '", coverage_file, "': ", e$message), call. = FALSE)
})
message("Coverage data read and processed successfully.")

# --- Calculate Weighted Chromosome Coverage ---
message("Calculating chromosome coverage...")
# Filter coverage data using the CHR list
chromosome_data <- coverage_data %>%
  dplyr::filter(contig_name %in% chr_names_list)

if (nrow(chromosome_data) == 0) {
    warning("Warning: No chromosome contigs from the list found in coverage data. Chromosome average depth will be 0.", call. = FALSE)
    chromosome_avg_depth <- 0
    chromosome_summary <- list(contig_count = 0, total_length = 0) # Dummy summary
} else {
    chromosome_summary <- chromosome_data %>%
      dplyr::summarise(
        total_length = sum(length, na.rm = TRUE),
        total_weighted_depth = sum(meandepth * length, na.rm = TRUE),
        contig_count = n(),
        .groups = 'drop'
      )

    if (chromosome_summary$total_length <= 0) {
      warning("Total chromosome length calculated is zero or negative. Chromosome average depth will be 0.", call. = FALSE)
      chromosome_avg_depth <- 0
    } else {
      chromosome_avg_depth <- chromosome_summary$total_weighted_depth / chromosome_summary$total_length
    }
}
message(sprintf("Chromosome Summary: %d contigs processed, Total Length %d, Weighted Avg Depth %.2f",
                chromosome_summary$contig_count,
                chromosome_summary$total_length,
                chromosome_avg_depth))


# --- Calculate Copy Number for Each Plasmid ---
message("Calculating plasmid copy numbers...")
# Filter coverage data to find plasmid rows
all_plasmid_coverage_data <- coverage_data %>%
  dplyr::filter(contig_name %in% plasmid_names_list)

results_list <- list()

if (length(plasmid_names_list) == 0) {
    warning("Plasmid name list is empty. No copy numbers to calculate.", call. = FALSE)
} else {
    if (nrow(all_plasmid_coverage_data) == 0) {
         warning("Warning: None of the plasmid names found in the list were present in the coverage summary file.", call. = FALSE)
    }
    # Iterate over each plasmid name in the list
    for (p_name in plasmid_names_list) {
      plasmid_row <- all_plasmid_coverage_data %>% dplyr::filter(contig_name == p_name)

      if (nrow(plasmid_row) == 1) {
        p_depth <- plasmid_row$meandepth[1]
        p_length <- plasmid_row$length[1]
        copy_number <- if (chromosome_avg_depth > 0) {
          p_depth / chromosome_avg_depth
        } else {
          NA # Assign NA if chromosome depth is zero
        }
        results_list[[p_name]] <- data.frame(
          plasmid_contig = p_name,
          length = p_length,
          mean_depth = p_depth,
          chromosome_mean_depth = chromosome_avg_depth,
          estimated_copy_number = copy_number
        )
      } else if (nrow(plasmid_row) == 0) {
          # Plasmid from list not found in coverage data
          warning(paste("Plasmid contig '", p_name, "' from list not found in coverage summary. Skipping."), call. = FALSE)
          results_list[[p_name]] <- data.frame(
            plasmid_contig = p_name, length = NA, mean_depth = NA,
            chromosome_mean_depth = chromosome_avg_depth, estimated_copy_number = NA
          )
      } else {
          # Multiple rows for the same contig name
          warning(paste("Multiple coverage entries found for plasmid contig '", p_name, "'. Using first entry."), call. = FALSE)
          p_depth <- plasmid_row$meandepth[1]
          p_length <- plasmid_row$length[1]
          copy_number <- if (chromosome_avg_depth > 0) p_depth / chromosome_avg_depth else NA
          results_list[[p_name]] <- data.frame(
              plasmid_contig = p_name, length = p_length, mean_depth = p_depth,
              chromosome_mean_depth = chromosome_avg_depth, estimated_copy_number = copy_number
          )
      }
    }
}

# --- Combine results and Write Output ---
if (length(results_list) > 0) {
  final_report <- dplyr::bind_rows(results_list)

  message(paste("Writing final report to:", output_file))
  tryCatch({
      readr::write_tsv(final_report, output_file)
      message("Output file written successfully.")
  }, error = function(e) {
      stop(paste("Error writing output file '", output_file, "': ", e$message), call. = FALSE)
  })

} else {
  warning("No plasmid copy numbers were calculated. Output file will be empty or may not be created.", call. = FALSE)
  # Create an empty file with headers
  tryCatch({
      # Define headers matching the data frame structure
      empty_df <- data.frame(plasmid_contig=character(), length=integer(), mean_depth=double(),
                             chromosome_mean_depth=double(), estimated_copy_number=double())
      readr::write_tsv(empty_df, output_file)
      message(paste("Empty report file with headers written to:", output_file))
  }, error = function(e){
      stop(paste("Error writing empty output file '", output_file, "': ", e$message), call. = FALSE)
  })
}

message("R script finished.")

