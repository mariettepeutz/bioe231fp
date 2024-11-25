#!/bin/bash

""" Instructions for running the script:

1. Ensure you have created the APACHE_ROOT environment variable pointing to your Apache root directory.
   Example: export APACHE_ROOT='/path/to/apache/root'
   If you're unsure of the correct path, you can find it by running:
   sudo find / -name "www" 2>/dev/null

2. Navigate to a temporary working directory (e.g., ~/tmp) to clone the repository.
   mkdir ~/tmp (if you haven't already)
   cd ~/tmp

3. Clone the repository containing the script:
   git clone https://github.com/mariettepeutz/bioe231fp.git

4. Move into the repository directory:
   cd bioe231fp

5. Make the script executable:
   chmod +x main.sh

6. Run the script:
   ./main.sh
   Ensure you have all required dependencies installed, such as wget, samtools, bowtie2, and jbrowse.

Additional Notes:
- The script assumes you have a JBrowse2 instance installed in the $APACHE_ROOT/jbrowse2 directory.
- If any errors occur, the script will exit and display an error message.
"""

"""Dengue"""

# Helper function for error handling - exits the code in case something goes wrong!
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error encountered in previous command. Exiting."
        exit 1
    fi
}

# Download and Process the Reference Genome #
echo "Downloading Dengue genome (reference genome)..."
# this is wrong!! wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" -O genome.fna.gz
check_error

echo "Decompressing reference genome file..."
gunzip genome.fna.gz
check_error

echo "Renaming reference genome file..."
mv genome.fna viral_genome.fa
check_error

echo "Indexing reference genome file with samtools..."
samtools faidx viral_genome.fa
check_error

echo "Building Bowtie2 index for reference genome..."
# Bowtie2 requires an indexed version of the reference genome for alignment
# The `bowtie2-build` command generates an index based on the reference genome
bowtie2-build viral_genome.fa viral_genome
check_error

echo "Adding reference genome assembly to JBrowse..."
jbrowse add-assembly viral_genome.fa --out /path/to/your/apache/root/jbrowse2 --load copy
check_error

# Download and Process Annotations #
echo "Downloading genome annotations..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz" -O annotations.gff.gz
check_error

echo "Decompressing annotations file..."
gunzip annotations.gff.gz
check_error

echo "Sorting GFF3 annotations..."
jbrowse sort-gff annotations.gff > genes.gff
check_error

echo "Compressing sorted GFF3 file..."
bgzip genes.gff
check_error

echo "Indexing compressed GFF3 file with tabix..."
tabix genes.gff.gz
check_error

echo "Adding annotations track to JBrowse..."
jbrowse add-track genes.gff.gz --out /path/to/your/apache/root/jbrowse2 --load copy
check_error

echo "Reference dengue genome and annotations successfully added to JBrowse."

# Align Comparison Genome to Reference Genome #
echo "Downloading comparison genome FASTA file..."
wget -q "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000862125.1/download?include_annotation_type=GENOME_FASTA&hydrated=FULLY_HYDRATED" -O comparison_genome.fna
check_error

echo "Aligning comparison genome to reference genome using Bowtie2..."
# Bowtie2 aligns the comparison genome (FASTA format) to the reference genome.
# The `-x viral_genome` specifies the reference genome index created earlier.
# The `-f comparison_genome.fna` specifies the input file (FASTA format).
# The output is saved as a SAM file (`comparison_genome.sam`), which contains the alignment data.
bowtie2 -x viral_genome -f comparison_genome.fna -S comparison_genome.sam
check_error

echo "Converting SAM file to BAM file..."
# Convert the SAM file (text-based) to BAM format (binary, optimized for storage and speed).
samtools view -bS comparison_genome.sam > comparison_genome.bam
check_error

echo "Sorting BAM file..."
# Sort the BAM file for proper visualization in JBrowse.
samtools sort comparison_genome.bam -o comparison_genome.sorted.bam
check_error

echo "Indexing BAM file..."
# Index the BAM file for quick access to specific regions during visualization.
samtools index comparison_genome.sorted.bam
check_error

echo "Adding alignment track to JBrowse..."
# Add the BAM file (alignment data) as a new track to JBrowse.
jbrowse add-track comparison_genome.sorted.bam --out /path/to/your/apache/root/jbrowse2 --load copy
check_error

echo "Comparison genome alignment successfully added to JBrowse."


#basic code
# # Upload main Dengue genome
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz

# gunzip GCF_000862125.1_ViralProj15306_genomic.fna.gz
# mv GCF_000862125.1_ViralProj15306_genomic.fna viral_genome.fa
# samtools faidx viral_genome.fa
# jbrowse add-assembly viral_genome.fa --out $APACHE_ROOT/jbrowse2 --load copy

# #Uploading annotations
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz
# gunzip GCF_000862125.1_ViralProj15306_genomic.gff.gz
# jbrowse sort-gff GCF_000862125.1_ViralProj15306_genomic.gff > genes.gff
# bgzip genes.gff
# tabix genes.gff.gz
# jbrowse add-track genes.gff.gz --out $APACHE_ROOT/jbrowse2 --load copy
