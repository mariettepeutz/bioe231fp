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

### Dengue ###

#!/bin/bash

# Helper function for error handling
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error encountered in previous command. Exiting."
        exit 1
    fi
}

### Define Directories ###
# Set the working directory to the bioe231fp folder
WORKDIR=$(pwd)  # Assumes the script is run from the bioe231fp folder
APACHE_ROOT=/opt/homebrew/var/www

echo "Working directory set to: $WORKDIR"

# Ensure Apache root is defined
if [ ! -d "$APACHE_ROOT/jbrowse2" ]; then
    echo "Error: JBrowse2 not found at $APACHE_ROOT/jbrowse2. Please check your APACHE_ROOT variable."
    exit 1
fi

### Download and Process the Reference Genome ###
echo "Downloading Dengue genome (reference genome)..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" -O $WORKDIR/genome.fna.gz
check_error

echo "Decompressing reference genome file..."
gunzip -f $WORKDIR/genome.fna.gz
check_error

echo "Renaming reference genome file..."
mv $WORKDIR/genome.fna $WORKDIR/viral_genome.fa
check_error

echo "Indexing reference genome file with samtools..."
samtools faidx $WORKDIR/viral_genome.fa
check_error

echo "Building Bowtie2 index for reference genome..."
bowtie2-build $WORKDIR/viral_genome.fa $WORKDIR/viral_genome
check_error

echo "Adding reference genome assembly to JBrowse..."
jbrowse add-assembly $WORKDIR/viral_genome.fa --out $APACHE_ROOT/jbrowse2 --load copy
check_error

### Download and Process Annotations ###
echo "Downloading genome annotations..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz" -O $WORKDIR/annotations.gff.gz
check_error

echo "Decompressing annotations file..."
gunzip -f $WORKDIR/annotations.gff.gz
check_error

echo "Sorting GFF3 annotations..."
jbrowse sort-gff $WORKDIR/annotations.gff > $WORKDIR/genes.gff
check_error

echo "Compressing sorted GFF3 file..."
bgzip -f $WORKDIR/genes.gff
check_error

echo "Indexing compressed GFF3 file with tabix..."
tabix $WORKDIR/genes.gff.gz
check_error

echo "Adding annotations track to JBrowse..."
jbrowse add-track $WORKDIR/genes.gff.gz --out $APACHE_ROOT/jbrowse2 --load copy
check_error

echo "Reference dengue genome and annotations successfully added to JBrowse."

### Align Comparison Genome to Reference Genome ###
echo "Downloading comparison genome FASTA file..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" -O $WORKDIR/comparison_genome.fna
check_error

echo "Aligning comparison genome to reference genome using Bowtie2..."
bowtie2 -x $WORKDIR/viral_genome -f $WORKDIR/comparison_genome.fna -S $WORKDIR/comparison_genome.sam
check_error

echo "Converting SAM file to BAM file..."
samtools view -bS $WORKDIR/comparison_genome.sam > $WORKDIR/comparison_genome.bam
check_error

echo "Sorting BAM file..."
samtools sort $WORKDIR/comparison_genome.bam -o $WORKDIR/comparison_genome.sorted.bam
check_error

echo "Indexing BAM file..."
samtools index $WORKDIR/comparison_genome.sorted.bam
check_error

echo "Adding alignment track to JBrowse..."
jbrowse add-track $WORKDIR/comparison_genome.sorted.bam --out $APACHE_ROOT/jbrowse2 --load copy
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
