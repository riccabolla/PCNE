# Plasmid Copy Estimator
A simple tool to estimate the copy number of plasmid contigs relative to chromosomal contigs using raw sequencing reads and FASTA files derived from a Platon analysis.
## Introduction
This tool automates the process of:<br>
1) Creating a combined reference from Platon-identified chromosome and plasmid contigs.<br>
2) Indexing the reference using BWA.<br>
3) Aligning paired-end sequencing reads to the combined reference using BWA-MEM.<br>
4) Sorting and indexing the resulting alignment (BAM file) using Samtools.<br>
5) Calculating per-contig coverage using samtools coverage.Calculating the length-weighted average depth of chromosome contigs.<br>
6) Calculating the estimated copy number for each plasmid contig relative to the average chromosome depth using an R script.<br>
## Installation<br>
The recommended way to install Plasmid Copy Estimator is via Conda, preferably through the Bioconda channel (once the recipe is submitted and accepted).<br>
1) **Set up Conda Channels (if not already done):**<br>
You need bioconda and conda-forge channels, with conda-forge having higher priority.<br>
conda config --add channels defaults<br>
conda config --add channels bioconda<br>
conda config --add channels conda-forge<br>
conda config --set channel_priority strict<br>
2) **Create a new environment and install:**<br>
#Replace 'plasmid_cn_env' with your desired environment name<br>
conda create -n plasmid_cn_env plasmid-copy-estimator<br>
conda activate plasmid_cn_env<br>
(Note: Until the package is accepted into Bioconda, you would install from a local build using conda install --use-local plasmid-copy-estimator after building it yourself with conda build conda/)<br>
## Dependencies. <br>
The tool relies on the following software, which will be installed automatically by Conda:<br>
BWA (>=0.7.17 recommended)<br>
Samtools (>=1.9 recommended)<br>
R (>=4.0 recommended)<br>
R Packages: readr, dplyr<br>
## Usage<br>
plasmid-copy-estimator -c <chromosome.fasta> -p <plasmid.fasta> -r1 <reads_R1.fastq.gz> -r2 <reads_R2.fastq.gz> [-t <threads>] [-o <output_prefix>]

Arguments: <br>
* -c FILE: Path to chromosome contigs FASTA file (from Platon) [Required] <br>
* -p FILE: Path to plasmid contigs FASTA file (from Platon) [Required] <br>
* -r1 FILE: Path to forward reads (FASTQ format, can be gzipped) [Required] <br>
* -r2 FILE: Path to reverse reads (FASTQ format, can be gzipped) [Required] <br>
* -t INT: Number of threads to use for alignment and sorting (default: 1) [Optional] <br>
* -o STR: Prefix for output files (default: plasmid_cn_result) [Optional] <br>
* -h: Display help message (implicit) Example# Activate the conda environment first conda activate plasmid_cn_env

# Run the tool
plasmid-copy-estimator \ <br>
  -c platon_output/my_sample.chromosome.fasta \ <br>
  -p platon_output/my_sample.plasmid.fasta \ <br>
  -r1 input_reads/my_sample_R1.fastq.gz \ <br>
  -r2 input_reads/my_sample_R2.fastq.gz \ <br>
  -t 8 \ <br>
  -o my_sample_copy_num

OutputThe tool generates several intermediate files (reference, index files, BAM file, coverage summary). <br>
The main output is a TSV (Tab-Separated Values) file named **<output_prefix>_plasmid_copy_numbers.tsv**

Example output.tsv: <br>
| plasmid_contig |length | mean_depth |chromosome_mean_depth |estimated_copy_number |
| ---------------|-------|------------|----------------------|----------------------|
|plasmid_contig_ 1|54321 |152.75|31.45|4.86|
|plasmid_contig_2_IncFIB|9876|28.50|31.45|0.91|
|...|...|...|...|...| 

Columns: <br>
* **plasmid_contig**: Name of the plasmid contig (from the input plasmid FASTA).<br>
* **length**: Length of the plasmid contig in base pairs.<br>
* **mean_depth**: Average sequencing depth calculated for this plasmid contig.<br>
* **chromosome_mean_depth**: The single, length-weighted average depth calculated across all chromosomal contigs.<br>
* **estimated_copy_number**: The calculated copy number (mean_depth / chromosome_mean_depth).<br>

**License**<br>
This project is licensed under the MIT License - see the LICENSE file for details.<br>
**Contact** <br>
riccardo.bollini@hunimed.ey <br>
**Issues**<br>
Please report any issues or suggestions via the GitHub Issues page.<br>