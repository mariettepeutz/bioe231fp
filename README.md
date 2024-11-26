# The Flavi Fighters
## BioE231 Final Project by Smrithi Surender, Marriëtte Peutz, & Kate Barouch

This project provides a customizable database installer for JBrowse2, focused on the Flavivirus family, to support pandemic researchers by integrating genome alignments, annotations, and bioinformatics tools for rapid pathogen identification and evolutionary analysis.

## **Table of Contents**
1. [Instructions for Running the Script](#instructions-for-running-the-script)
2. [The Flavivirus Family](#the-flavivirus-family)
3. [Data Sources](#data-sources)
4. [Usage Guide](#usage-guide)
5. [Collaborator Contributions](#collaborator-contributions)

---

## Instructions for running the script:

1. Install all necessary dependencies:
   Run the following command to install `wget`, `samtools`, `bowtie2`, `jbrowse`, and other essential tools:
      bash requirements.txt

2. Set the APACHE_ROOT environment variable:
   Ensure you have created the `APACHE_ROOT` environment variable pointing to your Apache root directory.
   Example:
       export APACHE_ROOT='/path/to/apache/root'

   If you're unsure of the correct path, run the following command to locate your Apache root directory:
       sudo find / -name "www" 2>/dev/null

3. Create a temporary working directory (if you haven't already) and navigate to it:
      mkdir ~/tmp
      cd ~/tmp

4. Download and copy over JBrowse 2 into the apache2 root dir, setting the owner to the current user with chown:
      jbrowse create output_folder
      sudo mv output_folder $APACHE_ROOT/jbrowse2
      sudo chown -R $(whoami) $APACHE_ROOT/jbrowse2

5. Test your JBrowse2 installation by opening the browser and typing in:
   
      http://localhost/jbrowse2/

7. Clone the repository containing the script:
   
   git clone https://github.com/mariettepeutz/bioe231fp.git

9. Move into the repository directory:
   cd bioe231fp

10. Make the script executable:
   chmod +x main.sh

11. Run the script:
   ./main.sh

---

## Flavivirus Family

We chose to focus on the Flavivirus family due to the increasing global health threats posed by this virus family with climate warming expanding mosquito habitats and contributing to outbreaks in new regions.

### Dengue Virus (DENV): 
A major global health threat, with over 6,800 cases reported in the United States as of October 31, 2024—more than double the cases reported in 2023—driven by climate change and increasing prevalence in the Americas and Asia.

### Zika Virus (ZIKV): 
Recognized for causing severe congenital outcomes like microcephaly, it is included on the WHO priority list for diseases with epidemic potential, and there is currently no vaccine available.

### West Nile Virus (WNV): 
A widespread virus transmitted by Culex mosquitoes, with severe neurological complications in some cases, and no available vaccine, underscoring the importance of further research.

---

## **Data Sources**

### **External Data Sources**:
- Genome sequences and annotations were retrieved from the NCBI FTP repository.
- Specific data includes genomic FASTA files and GFF annotation files for Dengue, Zika, and West Nile viruses.

#### **Example Links**:
- **Dengue**: [NCBI Genome Reference](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/)

### **Self-Generated Data**:
- Alignments between reference and comparison genomes.
- Custom annotations generated using JBrowse2 tools.

---

## **Usage Guide**

### **Navigating JBrowse2**
- **Reading Frames**: Explore six reading frames (three forward and three reverse).
- **Annotation Tracks**: View genome annotations (e.g., genes, regulatory regions).
- **Alignment Tracks**: Analyze alignment tracks for insights into different viral strains.

### **Customizing Data Views**


---

## **Collaborator Contributions**

### **Smrithi Surender**:
- 

### **Marriëtte Peutz**:
- 

### **Kate Barouch**:
- 
