#!/bin/bash

### Dengue Alignments Tracks ###

# Helper function for error handling
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error encountered in previous command. Exiting."
        exit 1
    fi
}

### Check Directory and APACHE_ROOT ###
WORKDIR=$(pwd)

# Check if APACHE_ROOT is defined
if [ -z "$APACHE_ROOT" ]; then
    echo "Error: APACHE_ROOT is not defined. Please set it using the following command:"
    echo "export APACHE_ROOT='/path/to/apache/root'"
    echo "Refer to the instructions in the script header or requirements.txt for details."
    exit 1
fi

echo "Working directory set to: $WORKDIR"
echo "Using Apache root directory: $APACHE_ROOT"

# Ensure JBrowse2 directory exists
if [ ! -d "$APACHE_ROOT/jbrowse2" ]; then
    echo "Error: JBrowse2 not found in $APACHE_ROOT/jbrowse2."
    echo "Please ensure JBrowse2 is installed and configured correctly."
    exit 1
fi

# Function to process a reference genome
process_reference_genome() {
    virus_name="$1"
    genome_url="$2"
    echo "Processing $virus_name reference genome..."
    wget -q "$genome_url" -O "$WORKDIR/${virus_name}_genome.fna.gz"
    check_error

    echo "Unzipping $virus_name reference genome..."
    gunzip -f "$WORKDIR/${virus_name}_genome.fna.gz"
    check_error

    echo "Indexing $virus_name reference genome..."
    samtools faidx "$WORKDIR/${virus_name}_genome.fna"
    check_error

    echo "Adding $virus_name reference genome to JBrowse..."
    jbrowse add-assembly "$WORKDIR/${virus_name}_genome.fna" --out "$APACHE_ROOT/jbrowse2" --load copy
    check_error
}

# Function to process a reference annotation
process_reference_annotation() {
    virus_name="$1"
    annotation_url="$2"
    echo "Processing $virus_name annotations..."
    wget -q "$annotation_url" -O "$WORKDIR/${virus_name}_annotations.gff.gz"
    check_error

    echo "Unzipping $virus_name annotations..."
    gunzip -f "$WORKDIR/${virus_name}_annotations.gff.gz"
    check_error

    echo "Sorting $virus_name annotations..."
    jbrowse sort-gff "$WORKDIR/${virus_name}_annotations.gff" > "$WORKDIR/${virus_name}_genes.gff"
    check_error

    echo "Compressing $virus_name annotations..."
    bgzip -f "$WORKDIR/${virus_name}_genes.gff"
    check_error

    echo "Indexing $virus_name annotations..."
    tabix "$WORKDIR/${virus_name}_genes.gff.gz"
    check_error

    echo "Adding $virus_name annotations to JBrowse..."
    jbrowse add-track "$WORKDIR/${virus_name}_genes.gff.gz" --out "$APACHE_ROOT/jbrowse2" --load copy
    check_error
}

# Function to process comparison genome for alignment track
process_comparison_genome() {
    reference_virus_name="$1"
    comparison_virus_name="$2"
    comparison_genome_url="$3"

    echo "Building Bowtie2 index for $reference_virus_name genome..."
    bowtie2-build "$WORKDIR/${reference_virus_name}_genome.fna" "$WORKDIR/${reference_virus_name}_genome_index"
    check_error

    echo "Downloading $comparison_virus_name..."
    wget -q "$comparison_genome_url" -O "$WORKDIR/${comparison_virus_name}.fa"
    check_error

    echo "Aligning $comparison_virus_name to $reference_virus_name reference genome..."
    bowtie2 -x "$WORKDIR/${reference_virus_name}_genome_index" -f "$WORKDIR/${comparison_virus_name}.fa" -S "$WORKDIR/${comparison_virus_name}.sam"
    check_error

    echo "Converting $comparison_virus_name alignment to BAM format..."
    samtools view -bS "$WORKDIR/${comparison_virus_name}.sam" > "$WORKDIR/${comparison_virus_name}.bam"
    check_error

    echo "Sorting $comparison_virus_name BAM file..."
    samtools sort "$WORKDIR/${comparison_virus_name}.bam" -o "$WORKDIR/${comparison_virus_name}.sorted.bam"
    check_error

    echo "Indexing $comparison_virus_name sorted BAM file..."
    samtools index "$WORKDIR/${comparison_virus_name}.sorted.bam"
    check_error

    echo "Adding $comparison_virus_name alignment track to JBrowse..."
    jbrowse add-track "$WORKDIR/${comparison_virus_name}.sorted.bam" --out "$APACHE_ROOT/jbrowse2" --load copy --name "$comparison_virus_name"
    check_error
}

### Main Script Execution ###

# Process DENV-1 Genomes
process_reference_genome "DENV-1" "https://ftp.ncbi..."
process_reference_annotation "DENV-1" "https://ftp.ncbi..."

process_comparison_genome "DENV-1" \
    "https://ftp.ncbi...." \
    "Variant1" "https://..."
process_comparison_genome "DENV-1" \
    "https://ftp.ncbi...." \
    "Variant2" "https://..."
process_comparison_genome "DENV-1" \
    "https://ftp.ncbi...." \
    "Variant3" "https://..."

# Process DENV-2 Genomes
process_reference_genome "DENV-2" "https://ftp.ncbi..."
process_reference_annotation "DENV-2" "https://ftp.ncbi..."

process_comparison_genome "DENV-2" \
    "https://ftp.ncbi...." \
    "Variant1" "https://..."
process_comparison_genome "DENV-2" \
    "https://ftp.ncbi...." \
    "Variant2" "https://..."
process_comparison_genome "DENV-2" \
    "https://ftp.ncbi...." \
    "Variant3" "https://..."

# Process DENV-3 Genomes
process_reference_genome "DENV-3" "https://ftp.ncbi..."
process_reference_annotation "DENV-3" "https://ftp.ncbi..."

process_comparison_genome "DENV-3" \
    "https://ftp.ncbi...." \
    "Variant1" "https://..."
process_comparison_genome "DENV-3" \
    "https://ftp.ncbi...." \
    "Variant2" "https://..."
process_comparison_genome "DENV-3" \
    "https://ftp.ncbi...." \
    "Variant3" "https://..."

# Process DENV-4 Genomes
process_reference_genome "DENV-4" "https://ftp.ncbi..."
process_reference_annotation "DENV-4" "https://ftp.ncbi..."

process_comparison_genome "DENV-4" \
    "https://ftp.ncbi...." \
    "Variant1" "https://..."
process_comparison_genome "DENV-4" \
    "https://ftp.ncbi...." \
    "Variant2" "https://..."
process_comparison_genome "DENV-4" \
    "https://ftp.ncbi...." \
    "Variant3" "https://..."

echo "All reference genomes, annotations, and comparison genomes have been processed and added to JBrowse."
