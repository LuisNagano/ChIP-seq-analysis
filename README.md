# Comprehensive ChIP-Seq Analysis Pipeline

Welcome to the **Comprehensive ChIP-Seq Analysis Pipeline** repository! This project demonstrates a complete workflow for ChIP-Seq data analysis, from raw data processing to peak calling and downstream analysis. It is designed to showcase best practices and can serve as a robust portfolio project to attract attention from leading companies like Illumina, Qiagen, Bayern, and Google.

## Table of Contents

1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Pipeline Overview](#pipeline-overview)
4. [Installation](#installation)
5. [Data Preparation](#data-preparation)
6. [Alignment with Bowtie2](#alignment-with-bowtie2)
7. [Quality Control](#quality-control)
8. [Peak Calling with MACS2](#peak-calling-with-macs2)
9. [Post-Processing with SAMtools](#post-processing-with-samtools)
10. [Visualization with IGV](#visualization-with-igv)
11. [Annotation with HOMER](#annotation-with-homer)
12. [Motif Analysis](#motif-analysis)
13. [IDR Analysis](#idr-analysis)
14. [Usage](#usage)
15. [Contributing](#contributing)
16. [License](#license)
17. [Additional Recommendations](#additional-recommendations)

## Introduction

Chromatin Immunoprecipitation followed by sequencing (ChIP-Seq) is a powerful technique to investigate protein-DNA interactions, histone modifications, and transcription factor binding sites across the genome. This pipeline provides a standardized approach to process ChIP-Seq data, ensuring reproducibility and scalability.

## Prerequisites

Before you begin, ensure you have the following software and tools installed:

- **Bowtie2**: For aligning sequencing reads to a reference genome.
- **SAMtools**: For handling SAM/BAM files.
- **MACS2**: For peak calling.
- **HOMER**: For peak annotation and motif analysis.
- **IGV (Integrative Genomics Viewer)**: For visualizing alignment and peaks.
- **IDR (Irreproducible Discovery Rate)**: For assessing the reproducibility of peaks across replicates.
- **Other Dependencies**:
  - Perl
  - Ghostscript
  - WebLogo
  - Blat
  - liftOver
  - bedGraphToBigWig

## Pipeline Overview

The pipeline consists of the following major steps:

1. **Installation of Tools**: Setting up necessary software.
2. **Data Preparation**: Downloading reference genomes and input files.
3. **Alignment**: Mapping reads to the reference genome using Bowtie2.
4. **Quality Control**: Assessing alignment quality.
5. **Peak Calling**: Identifying enriched regions using MACS2.
6. **Post-Processing**: Converting and sorting SAM/BAM files with SAMtools.
7. **Visualization**: Inspecting alignments and peaks with IGV.
8. **Annotation**: Annotating peaks with HOMER.
9. **Motif Analysis**: Discovering enriched motifs within peaks.
10. **IDR Analysis**: Evaluating peak reproducibility across biological replicates.

## Installation

### Bowtie2

```bash
# Download Bowtie2
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download -O bowtie2.zip

# Unzip
unzip bowtie2.zip -d ~/ChIP-Bioinfo/Programs/Bowtie2

# Add Bowtie2 to PATH
echo 'export PATH="$HOME/ChIP-Bioinfo/Programs/Bowtie2/bowtie2-2.4.5:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### SAMtools

```bash
# Clone SAMtools repository
git clone https://github.com/samtools/samtools.git
cd samtools
./configure
make
sudo make install

# Verify installation
samtools --version
```

### MACS2

```bash
# Install MACS2 via pip
pip install MACS2

# Verify installation
macs2 --version
```

### HOMER

```bash
# Download HOMER
wget https://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
```

### Other Dependencies

Ensure the following tools are installed and added to your PATH:

- **Perl**: Required for HOMER scripts.
- **Ghostscript**: For HOMER motif analysis.
- **WebLogo**: For visualizing sequence logos.
- **Blat**, **liftOver**, **bedGraphToBigWig**: For various genomic data conversions.

**Example installation for Ghostscript:**

```bash
# Download Ghostscript
wget https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs952/ghostscript-9.52.tar.gz
tar -zxvf ghostscript-9.52.tar.gz
cd ghostscript-9.52
./configure
make
sudo make install
```

**Add other tools to PATH as needed:**

```bash
export PATH="$HOME/ChIP-Bioinfo/Programs/WebLogo/:$PATH"
export PATH="$HOME/ChIP-Bioinfo/Programs/blat/:$PATH"
export PATH="$HOME/ChIP-Bioinfo/Programs/liftOver/:$PATH"
export PATH="$HOME/ChIP-Bioinfo/Programs/bedGraphToBigWig/:$PATH"
```

## Data Preparation

1. **Download Reference Genome**: Obtain the reference genome for your organism (e.g., mouse `mm10` or human `hg38`) from [Ensembl](https://www.ensembl.org/index.html) or [UCSC Genome Browser](https://genome.ucsc.edu/).

2. **Organize Directory Structure**:

    ```bash
    mkdir -p ~/ChIP-Bioinfo/{Programs,Data/{Reference,Raw,Aligned,Peaks,Annotations,Visualization,IDR,Motifs}}
    ```

3. **Place Reference Genome**:

    ```bash
    cp path_to_reference_genome.fna ~/ChIP-Bioinfo/Data/Reference/
    ```

4. **Build Bowtie2 Index**:

    ```bash
    bowtie2-build ~/ChIP-Bioinfo/Data/Reference/mm10.fna ~/ChIP-Bioinfo/Data/Reference/mm10
    ```

## Alignment with Bowtie2

Align your ChIP-Seq and control (input) FASTQ files to the reference genome.

### Command Structure

```bash
bowtie2 -q <input_fastq> -x <bowtie2_index> --time -S <output_sam>
```

### Example

```bash
# For Control Sample
bowtie2 -q ~/ChIP-Bioinfo/Data/Raw/INPUT/SRR836134.fastq -x ~/ChIP-Bioinfo/Data/Reference/mm10 --time -S ~/ChIP-Bioinfo/Data/Aligned/INPUT_aligned.sam

# For ChIP Sample
bowtie2 -q ~/ChIP-Bioinfo/Data/Raw/ChIP_COUPTFII/SRR836133.fastq -x ~/ChIP-Bioinfo/Data/Reference/mm10 --time -S ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.sam
```

## Quality Control

Assess the quality of your alignments by reviewing Bowtie2 alignment statistics:

```bash
# Example Output Interpretation
32708392 reads; of these:
  32708392 (100.00%) were unpaired; of these:
    378740 (1.16%) aligned 0 times
    23870993 (72.98%) aligned exactly 1 time
    8458659 (25.86%) aligned >1 times
98.84% overall alignment rate
```

High alignment rates (typically >90%) indicate good quality alignments.

## Peak Calling with MACS2

Identify enriched regions (peaks) using MACS2.

### Command Structure

```bash
macs2 callpeak -t <ChIP_SAM> -c <Control_SAM> -f SAM -g <genome_size> -p <pvalue_cutoff> -B -n <output_name>
```

### Example

```bash
macs2 callpeak -t ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.sam \
              -c ~/ChIP-Bioinfo/Data/Aligned/INPUT_aligned.sam \
              -f SAM \
              -g mm \
              -p 0.05 \
              -B \
              -n MACS2_COUPTFII
```

**Arguments Explained**:

- `-t`: ChIP-Seq aligned SAM file.
- `-c`: Control (input) aligned SAM file.
- `-f`: Format of input files (`SAM` in this case).
- `-g`: Genome size (`mm` for mouse, `hs` for human).
- `-p`: P-value cutoff for peak significance.
- `-B`: Generate bedGraph files for visualization.
- `-n`: Name prefix for output files.

## Post-Processing with SAMtools

Convert SAM files to BAM, sort, and index for downstream applications.

### Convert SAM to BAM

```bash
samtools view -bS ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.sam > ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.bam
```

### Sort BAM Files

```bash
samtools sort ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.bam -o ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.sorted.bam
samtools sort ~/ChIP-Bioinfo/Data/Aligned/INPUT_aligned.bam -o ~/ChIP-Bioinfo/Data/Aligned/INPUT_aligned.sorted.bam
```

### Index BAM Files

```bash
samtools index ~/ChIP-Bioinfo/Data/Aligned/ChIP_COUPTFII_aligned.sorted.bam
samtools index ~/ChIP-Bioinfo/Data/Aligned/INPUT_aligned.sorted.bam
```

## Visualization with IGV

Use [IGV](https://software.broadinstitute.org/software/igv/) to visualize alignments and peaks.

1. **Load BAM and BAI Files**: Ensure `.bam` and `.bam.bai` files are in the same directory.
2. **Load Peaks**: Import MACS2 peak files (`.narrowPeak` or `.bed`).

## Annotation with HOMER

Annotate peaks to identify associated genes and genomic features.

### Convert MACS2 Peaks to BED

Extract chromosome, start, end, and peak ID from MACS2 output:

```bash
awk '{print $1"\t"$2"\t"$3"\t"$4}' MACS2_COUPTFII_peaks.narrowPeak > COUPTFII_peaks.bed
```

### Install Reference Genome in HOMER

```bash
perl ~/ChIP-Bioinfo/Programs/HOMER/configureHomer.pl -install mm10
```

### Annotate Peaks

```bash
annotatePeaks.pl ~/ChIP-Bioinfo/Data/Peaks/COUPTFII_peaks.bed mm10 > ~/ChIP-Bioinfo/Data/Annotations/COUPTFII_annotation.txt -go ~/ChIP-Bioinfo/Data/Annotations/ -genomeOntology ~/ChIP-Bioinfo/Data/Annotations/
```

## Motif Analysis

Discover enriched DNA motifs within peaks using HOMER.

### Prepare BED File

Ensure the BED file contains chromosome, start, end, and peak ID.

### Run HOMER Motif Discovery

```bash
findMotifsGenome.pl ~/ChIP-Bioinfo/Data/Peaks/COUPTFII_peaks.bed mm10 ~/ChIP-Bioinfo/Data/Motifs/COUPTFII_motifs -size 200 -mask -oligo
```

**Arguments Explained**:

- `-size`: Window size around peaks (200 bp recommended for transcription factors).
- `-mask`: Masks low-complexity regions to reduce false positives.
- `-oligo`: Generates files with discovered motifs.

## IDR Analysis

Assess the reproducibility of peaks across biological replicates using IDR.

### Install IDR

```bash
# Clone IDR repository
git clone https://github.com/spendhirk/idr.git
cd idr
python setup.py install
```

### Run IDR

```bash
idr --samples ~/ChIP-Bioinfo/Data/Peaks/rep1_COUPTFII_peaks.narrowPeak ~/ChIP-Bioinfo/Data/Peaks/rep2_COUPTFII_peaks.narrowPeak \
    --input-file-type narrowPeak \
    --output-file ~/ChIP-Bioinfo/Data/IDR/COUPTFII_IDR_results.txt
```

## Usage

Clone the repository and follow the step-by-step instructions to execute the pipeline.

```bash
git clone https://github.com/yourusername/chipseq-analysis-pipeline.git
cd chipseq-analysis-pipeline
```

Refer to the [Installation](#installation) and [Pipeline Overview](#pipeline-overview) sections to set up and run the analysis.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.
