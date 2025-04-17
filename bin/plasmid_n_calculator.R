#!/usr/bin/env Rscript

# Script to calculate copy number for multiple plasmids relative to chromosome(s)
# Takes samtools coverage output and the original plasmid fasta file as input.

# --- Load Libraries ---
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
# Optional: library(Biostrings) or library(seqinr) for more robust FASTA parsing
# Using base R string manipulation here for fewer dependencies

# --- Get Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  message("Usage: Rscript calculate_multi_plasmid_cn.R <coverage_summary.txt> <plasmid.fasta> <output.tsv>")
  message("  <coverage_summary.txt>: Output file from 'samtools coverage'")
  message("  <plasmid.fasta>: Original plasmid FASTA file (used to identify plasmid contigs)")
  message("  <output.tsv>: Path for the final TSV report")
  stop("Incorrect number of arguments supplied.", call. = FALSE)
}

coverage_file <- args[1]
plasmid_fasta_file <- args[2]
output_file <- args[3]

# --- Function to Parse FASTA Headers ---
# Basic parser, assumes header is ">contig_name description" or just ">contig_name"
get_fasta_headers <- function(fasta_file) {
  tryCatch({
    lines <- readLines(fasta_file)
    header_lines <- lines[startsWith(lines, ">")]
    # Remove ">" and take only the first part before any space
    headers <- sub("^>", "", header_lines)
    headers <- sapply(strsplit(headers, "\\s+"), `[`, 1)
    if (length(headers) == 0) {
        stop("No headers found in plasmid FASTA file.", call. = FALSE)
    }
    return(unique(headers))
  }, error = function(e) {
    stop(paste("Error reading plasmid FASTA file '", fasta_file, "': ", e$message), call. = FALSE)
  })
}

# --- Read Plasmid Contig Names ---
message(paste("Reading plasmid contig names from:", plasmid_fasta_file))
plasmid_names_list <- get_fasta_headers(plasmid_fasta_file)
message(paste("Found", length(plasmid_names_list), "plasmid contig names."))
if (length(plasmid_names_list) == 0) stop("No plasmid contigs identified from FASTA.", call.=FALSE)

# --- Read and Process Coverage Data ---
message(paste("Reading coverage file:", coverage_file))
tryCatch({
  coverage_data <- readr::read_tsv(coverage_file, comment = "#", col_names = TRUE, show_col_types = FALSE) %>%
    {
      required_cols <- c("#rname", "startpos", "endpos", "meandepth")
      if (!all(required_cols %in% names(.))) {
        missing_cols <- setdiff(required_cols, names(.))
        stop(paste("Coverage file missing required columns:", paste(missing_cols, collapse=", ")), call. = FALSE)
      }
      .
    } %>%
    dplyr::rename(contig_name = `#rname`) %>%
    dplyr::mutate(length = endpos - startpos,
                  across(c(startpos, endpos, length, meandepth), as.numeric)) %>%
    dplyr::filter(length > 0) %>%
    { if(nrow(.) == 0) stop("No valid contigs found in coverage data after filtering.", call. = FALSE); . }

}, error = function(e) {
  stop(paste("Error reading/processing coverage file '", coverage_file, "': ", e$message), call. = FALSE)
})

# --- Separate Chromosome and Plasmid Data ---
# Chromosome data: contigs NOT in the plasmid list from the FASTA file
chromosome_data <- coverage_data %>%
  dplyr::filter(!contig_name %in% plasmid_names_list)

# Plasmid data: contigs that ARE in the plasmid list
all_plasmid_data <- coverage_data %>%
  dplyr::filter(contig_name %in% plasmid_names_list)

# --- Validate Data ---
if (nrow(chromosome_data) == 0) {
    warning("Warning: No chromosome contigs found in coverage data (all contigs matched plasmid names?). Check inputs.", call. = FALSE)
    # Allow proceeding, but copy number will be NA or Inf
    chromosome_avg_depth <- 0
    chromosome_summary <- list(contig_count = 0, total_length = 0) # Dummy summary
} else {
    # --- Calculate Weighted Chromosome Coverage ---
    chromosome_summary <- chromosome_data %>%
      dplyr::summarise(
        total_length = sum(length, na.rm = TRUE),
        total_weighted_depth = sum(meandepth * length, na.rm = TRUE),
        contig_count = n(),
        .groups = 'drop'
      )

    if (chromosome_summary$total_length <= 0) {
      warning("Total chromosome length is zero. Cannot calculate chromosome average depth.", call. = FALSE)
      chromosome_avg_depth <- 0
    } else {
      chromosome_avg_depth <- chromosome_summary$total_weighted_depth / chromosome_summary$total_length
    }
}

message(sprintf("Chromosome Summary: %d contigs, Total Length %d, Weighted Avg Depth %.2f",
                chromosome_summary$contig_count,
                chromosome_summary$total_length,
                chromosome_avg_depth))

# --- Calculate Copy Number for Each Plasmid ---
results_list <- list()

if (nrow(all_plasmid_data) == 0 && length(plasmid_names_list) > 0) {
    warning("Warning: None of the plasmid names found in FASTA were present in the coverage summary file.", call. = FALSE)
}

for (p_name in plasmid_names_list) {
  plasmid_row <- all_plasmid_data %>% dplyr::filter(contig_name == p_name)

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
      warning(paste("Plasmid contig '", p_name, "' from FASTA not found in coverage summary. Skipping."), call. = FALSE)
      # Optionally add a row indicating it was missing
      results_list[[p_name]] <- data.frame(
        plasmid_contig = p_name,
        length = NA, mean_depth = NA, chromosome_mean_depth = chromosome_avg_depth, estimated_copy_number = NA
      )
  } else {
      warning(paste("Multiple entries found for plasmid contig '", p_name, "' in coverage summary? Using first."), call. = FALSE)
      # Handle multiple entries if necessary, here just taking the first
      p_depth <- plasmid_row$meandepth[1]
      p_length <- plasmid_row$length[1]
      copy_number <- if (chromosome_avg_depth > 0) p_depth / chromosome_avg_depth else NA
      results_list[[p_name]] <- data.frame(
          plasmid_contig = p_name, length = p_length, mean_depth = p_depth,
          chromosome_mean_depth = chromosome_avg_depth, estimated_copy_number = copy_number
      )
  }
}

# --- Combine results and Write Output ---
if (length(results_list) > 0) {
  final_report <- dplyr::bind_rows(results_list)

  message(paste("Writing final report to:", output_file))
  readr::write_tsv(final_report, output_file)
} else {
  warning("No plasmid copy numbers were calculated. Output file will be empty or not created.", call. = FALSE)
  # Create an empty file or a file with just headers
  file.create(output_file) # Or write headers: write_tsv(data.frame(...), output_file)
}

message("R script finished.")

