### Dengue JBrowse Database Code ###

# Helper function for error handling - exits the code in case something goes wrong!
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error encountered in previous command. Exiting."
        exit 1
    fi
}

# Download and process the main genome (DENV-2 with V91A mutation)
echo "Downloading Dengue genome..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" -O genome.fna.gz
check_error

echo "Decompressing genome file..."
gunzip genome.fna.gz
check_error

echo "Renaming genome file..."
mv genome.fna viral_genome.fa
check_error

echo "Indexing genome file with samtools..."
samtools faidx viral_genome.fa
check_error

echo "Adding genome assembly to JBrowse..."
jbrowse add-assembly viral_genome.fa --out /path/to/your/apache/root/jbrowse2 --load copy
check_error

# Download and process annotations
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

echo "Genome and annotations successfully added to JBrowse for main Dengue genome (DENV-2 with V91A mutation)."

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
