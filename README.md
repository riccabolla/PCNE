![Static Badge](https://img.shields.io/badge/License-MIT-blue)
![Static Badge](https://img.shields.io/badge/Version-0.1.0-blue)

# Plasmid Copy Number Estimator
**P**lasmid **C**opy **N**umber **E**stimator (**PCNE**) is a simple tool to estimate the copy number of plasmid from an assembled genome.
## Disclaimer
This tool requires a previous step of plasmid identification using tools like Platon (recommended), MOB-Suite, PlasmidFinder... <br>
## Introduction
This tool automates the process of:<br>
1) Creating a combined reference from chromosome and plasmid contigs.<br>
2) Indexing the reference using BWA.<br>
3) Aligning paired-end sequencing reads to the combined reference using BWA-MEM.<br>
4) Sorting and indexing the resulting alignment (BAM file) using Samtools.<br>
5) Calculating per-contig coverage using samtools coverage. <br>
6) Calculating the length-weighted average depth of chromosome contigs.<br>
7) Calculating the estimated copy number for each plasmid contig relative to the average chromosome depth using an R script.<br>
## Installation<br>
[![Static Badge](https://img.shields.io/badge/Install_with-Bioconda-blue)](https://bioconda.github.io/)
The recommended way to install Plasmid Copy Number Estimator is via [BioConda](https://bioconda.github.io/).<br>
1) **Set up Conda Channels (if not already done):**<br>

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
2) **Create a new environment and install:**<br>
```
conda create -n pcne_env -c conda-forge -c bioconda pcne
conda activate pcne_env
```
## Dependencies. <br>
The tool relies on the following softwares, which will be installed automatically by Conda:<br>
1) **BWA** (>=0.7.18 recommended)<br>
2) **Samtools** (>=1.2 recommended)<br>
3) **R** (>=4.4.2 recommended)<br>
4) **R Packages**: readr (>=2.1.5), dplyr (>=1.1.4)<br>
## Usage<br>
```
pcne -c <chromosome.fasta> -p <plasmid.fasta> -r <reads_R1.fastq.gz> -R <reads_R2.fastq.gz> [-t <threads>] [-o <output_prefix>]
```
## Command line options: <br>
```
-c FILE       Path to chromosome contigs FASTA file (mode 1)
-p FILE       Path to plasmid contigs FASTA file (mode 1)
-a FILE       Path to assembled genome (FASTA file) (mode 2)
-C FILE       Path to the chromosome contigs list (mode 2)
-P FILE       Path to the plasmid contigs list (mode 2)
-r FILE       Path to forward reads (FASTQ format, can be gzipped) [Required] 
-R FILE       Path to reverse reads (FASTQ format, can be gzipped) [Required] 
-t INT        Number of threads to use for alignment and sorting (default: 1) [Optional] 
-o STR        Prefix for output files (default: pcne_result) [Optional] 
-h            Display help message (implicit) 
```

# Run the tool

The tool can be run in two different ways: <br>
**Mode 1**: it requires two separate `FASTA` files for chromosome and plasmid(s). <br>
```
#Example Mode 1
pcne \ 
  -c my_sample.chromosome.fasta \ 
  -p my_sample.plasmid.fasta \ 
  -r my_sample_R1.fastq.gz \ 
  -R my_sample_R2.fastq.gz \ 
  -t 8 \ 
  -o my_sample_pcne
```
**Mode 2**: it requires an assembled `FASTA` file, a list file with contig(s) assigned to chromosome, and a list file of contig(s) assigned to plasmid(S).
The list should be structured as follow:
```
plasmid1_contig
plasmid2_contig
plasmid3_contig
...
```
```
#Example Mode 2
pcne \ 
  -a my_sample_assembly.fasta
  -C chromosome.list
  -P plasmid.list \ 
  -r my_sample_R1.fastq.gz \ 
  -R my_sample_R2.fastq.gz \ 
  -t 8 \ 
  -o my_sample_pcne
```
**Note**: if files are not in the working folder, provide the PATH. <br>

The tool generates several intermediate files (reference, index files, BAM file, coverage summary). <br>
For both modes the main output is a `TSV` file.
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
riccardo.bollini@hunimed.eu <br>
**Issues**<br>
Please report any issues or suggestions via the GitHub Issues page.<br>
