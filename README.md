![Static Badge](https://img.shields.io/badge/License-MIT-blue)
![Static Badge](https://img.shields.io/badge/version-0.2.0-blue)

# Plasmid Copy Number Estimator
**P**lasmid **C**opy **N**umber **E**stimator (**PCNE**) is a simple tool to estimate the copy number of plasmid from an assembled genome. <br>
## ⚠️ Warning⚠️ 
This tool is currently in its early development stage (v0.2.0). While functional for its core purpose, it will receive updates often, including additional features, bug fixing... (for more information see [here](#Next-features)). Please report any issues or suggestions!
## Introduction
Determining the copy number of plasmids relative to the host chromosome is essential for understanding plasmid biology, evolution, and the dosage of plasmid-borne genes (e.g., antimicrobial resistance genes). PCNE automates this estimation from standard bioinformatics file formats. <br>
### Key features
* **Flexible input**: Accepts either pre-separated chromosome and plasmid FASTA files or a complete genome assembly FASTA with corresponding contig lists. It also allows to use as plasmid input a multi-fasta file with one contig per plasmid, a complete assembled plasmid (1 contig), or a draft assembled plasmid (one plasmid with multiple contigs). 
* **Multiple normalization**: To assess baseline coverage depth, it allows for a standard length-weighted average coverage of all designated chromosomal contigs or a more precise normalization using single-copy core genes identified bu BUSCO. 
* **Alignment Filtering**: Allow to filter alignments based on mapping quality and SAM flags before coverage calculation, helping to remove ambiguous or low-quality mappings.
* **Visualization**: Optionally generates a bar plot visualizing the estimated copy numbers for each plasmid

## Installation
### Bioconda [![Static Badge](https://img.shields.io/badge/Install_with-Bioconda-blue)](https://bioconda.github.io/)

The recommended way to install Plasmid Copy Number Estimator will be via [BioConda](https://bioconda.github.io/) once approved.<br>
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
### Source
Clone this repository to install the latest version direct from GitHub
```
cd $HOME
git clone https://github.com/riccabolla/PCNE.git 
cd PCNE/
```
(Note: While the Bioconda submission is pending, you can install it from your local build using `conda install --use-local pcne` after building it yourself with `conda build conda/`)
## Dependencies
The tool relies on the following softwares, which will be installed automatically by Conda:<br>
1) **BWA** (>=0.7.18 recommended)<br>
2) **Samtools** (>=1.2 recommended)<br>
3) **Busco** (=5.8.2)
4) **bedtools** (>=2.31.1)
5) **R** (>=4.4.2 recommended)<br>
6) **R Packages**: readr (>=2.1.5), dplyr (>=1.1.4), ggplot2(>=3.5.1)<br>
## Requirements
This tool requires a previous step of plasmid identification using tools like Platon (recommended), MOB-Suite, PlasmidFinder...
## Usage
```
pcne -c <chromosome.fasta> -p <plasmid.fasta> -r <reads_R1.fastq.gz> -R <reads_R2.fastq.gz> [-t <threads>] [-o <output_prefix>]
```
## Command line options
```
  -c, --chromosome FILE      Path to chromosome FASTA file (Required)  
  -p, --plasmid FILE         Path to plasmid FASTA file (Required)  
                             Use with `--single-plasmid` if file contains one fragmented plasmid  
  -a, --assembly FILE        Path to the assembled genome FASTA file (Required)  
  -C, --chr-list FILE        Path to file containing chromosome contig names (Required)  
  -P, --plasmid-list FILE    Path to file containing plasmid contig names (Required)  
  -r, --reads1 FILE          Path to forward reads (FASTQ) (Mandatory)  
  -R, --reads2 FILE          Path to reverse reads (FASTQ) (Mandatory)  
  -b, --busco                Activate BUSCO SCG normalization (default: OFF)  
  -L, --busco-lineage STR    BUSCO lineage dataset name (default: bacteria_odb12)  
  -Q, --min-quality INT      Minimum mapping quality (MQ) for read filtering (default: OFF)  
  -F, --filter INT           SAM flag to exclude reads (default: OFF)  
  -l, --plot                 Generate a plot of estimated copy numbers (.png)  
  -s, --single-plasmid       Treat all contigs in `-p` FASTA as one fragmented plasmid (Mode 1 only)  
  -t, --threads INT          Number of threads to use (default: 1)  
  -o, --output STR           Prefix for output files (default: pcne)  
  -k, --keep-intermediate    Keep intermediate files (default: OFF)  
  -v, --version              Show version information  
  -h, --help                 Show help message 
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
  -a my_sample_assembly.fasta \
  -C chromosome.list \
  -P plasmid.list \ 
  -r my_sample_R1.fastq.gz \ 
  -R my_sample_R2.fastq.gz \ 
  -t 8 \ 
  -o my_sample_pcne
```
**Note**: if files are not in the working folder, provide the PATH. <br>

For both modes the main output is a `TSV` file. <br>
Example `output.tsv`: <br>

| plasmid_contig |length | mean_depth |baseline_mean_depth |normalization_mode |estimated_copy_number |
|---|---|---|---|---|---|
|plasmid_contig_ 1|54321 |152.75|31.45|BUSCO_SCG|4.86|
|plasmid_contig_2_IncFIB|9876|28.50|31.45|BUSCO_SCG|0.91|
|...|...|...|...|...|...| 

Columns: <br>
* **plasmid_contig**: Name of the plasmid contig (from the input plasmid FASTA).<br>
* **length**: Length of the plasmid contig in base pairs.<br>
* **mean_depth**: Average sequencing depth calculated for this plasmid contig.<br>
* **baseline_mean_depth**: Baseline coverage depth.<br>
* **normalization mdoe**: how baseline coverage depth was calculated <br>
* **estimated_copy_number**: The calculated copy number (mean_depth / baseline_mean_depth).<br>

## <a name="Next-features"></a>Next features
### Major updates
* GC normalization

### **License**<br>
This project is licensed under the MIT License - see the [LICENSE](https://github.com/riccabolla/PCNE/blob/main/LICENSE) file for details.<br>

### **Contact** <br>
riccardo.bollini@hunimed.eu <br>

### **Issues**<br>
Please report any issues or suggestions via the GitHub [Issues](https://github.com/riccabolla/PCNE/issues) page.<br>
