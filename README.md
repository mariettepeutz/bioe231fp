# The Dengue Destroyers

This project provides a customizable database installer for JBrowse2, focused on the Dengue virus, to support pandemic researchers by integrating genome alignments and annotations for evolutionary analysis.

## **Table of Contents**
1. [Instructions for Running the Script](#instructions-for-running-the-script)
2. [Dengue Virus Background Information](#Dengue-Virus(DENV))
3. [Data Sources](#data-sources)
4. [JBrowse Usage Guide](#JBrowse-usage-guide)
5. [Collaborator Contributions](#collaborator-contributions)

---

## Instructions for running the script:

1. Clone the repository containing the script:
   
        git clone https://github.com/mariettepeutz/bioe231fp.git

2. Move into the repository directory:

        cd bioe231fp
   
4. Install all necessary dependencies:
   Based on your operating system, run the appropriate command to install wget, samtools, bowtie2, jbrowse, and other essential tools:

   For macOS:

       bash mac_requirements.sh
   For Linux:
   
       bash linux_requirements.sh

6. Set the APACHE_ROOT environment variable:
   Ensure you have created the `APACHE_ROOT` environment variable pointing to your Apache root directory.

      Example:
   
       export APACHE_ROOT='/path/to/apache/root'

   For a normal linux installation, the folder should be /var/www or /var/www/html, whereas when you install on macOS using brew it will likely be in /opt/homebrew/var/www (for M1) or /usr/local/var/www (for Intel). You can run brew --prefix to get the brew install location, and then from there it is in the var/www folder.

   If you're unsure of the correct path, run the following command to locate your Apache root directory:
   
       sudo find / -name "www" 2>/dev/null

7. Download and copy over JBrowse 2 into the apache2 root dir, setting the owner to the current user with chown:

       jbrowse create output_folder
       sudo mv output_folder $APACHE_ROOT/jbrowse2
       sudo chown -R $(whoami) $APACHE_ROOT/jbrowse2

8. Test your JBrowse2 installation by opening the browser and typing in:

      For macOS:
      
        http://localhost:8080/jbrowse2/
   
      For AWS: (xx.xxx.xxx.xx is your IP address)
      
        http://xx.xxx.xxx.xx/jbrowse2/
   
      For Linux: (xx.xxx.xxx.xx is your IP address)
      
        http://xx.xxx.xxx.xx:8080/jbrowse2/

10. Make the script executable:

        chmod +x main.sh

11. Run the script:

        ./main.sh

12. Navigate again to your JBrowse installation to see your final visualization of the viral genomes:

      For macOS:
      
        http://localhost:8080/jbrowse2/
   
      For AWS: (xx.xxx.xxx.xx is your IP address)
      
        http://xx.xxx.xxx.xx/jbrowse2/
   
      For Linux: (xx.xxx.xxx.xx is your IP address)
      
        http://xx.xxx.xxx.xx:8080/jbrowse2/

---

## Dengue Virus(DENV)

An increasing and major global health threat, with over 6,800 cases reported in the United States as of October 31, 2024—more than double the cases reported in 2023—driven by climate change expanding mosquito habitats and contributing to outbreaks in new regions, specifically in the Americas and Asia.

---


## Data Sources

### **External Data Sources**:
- **Genome sequences** (FASTA files) were retrieved from the NCBI FTP repository. For each sequence, we included NCBI Accession IDs, links and URLs to the published articles they were described in, as part of the Supplementary Data 1 .csv file in this repository.
- **Reference genomes and annotation files** (GFF files) were sourced from NCBI Genome. Documentation on this sourcing is in the Supplementary Data 2 PDF file in this repository. We manually edited to remove the two annotations for the polygene and the cds, as well as any references to them as 'parent' in the other annotations. Without this, JBrowse2 will only display the Polygene annotation. 
- Additional literature consulted is included in the Thematic Focus Write-up PDF in this repository. 

### **Self-Generated Data**:
- Alignments between reference and comparison genomes, using bowtie2.

---


### **Navigating JBrowse2**
- **Reading Frames**: Explore six reading frames (three forward and three reverse).
- **Annotation Tracks**: View genome annotations (e.g., genes, regulatory regions). Note the annotations for each protein are called mature_protein_regions. We have included an annotation for the V91A mutation in the DENV-3 strain, at position 7316.
- **Alignment Tracks**: Analyze alignment tracks for insights into different viral strains.

### **Customizing Data Views**
- Users can "Add" > "Linear Genome View" for each serotype, and add the alignment tracks they wish to visualize in the "Open Track Selector". Using the cursor, they can diminish space between alignment tracks.
- Users can align other sequences on the NCBI Genbank by calling the process_comparison_genome function in our main.sh on the serotype, a descriptive name, and the NCBI link directly to the fasta file. 
A command looks like this: 

      process_comparison_genome "DENV-4" \
          "DENV-4_Descriptive_Name" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=ACCESSION_ID&report=fasta&format=text" 


---

## Collaborator Contributions
- Smrithi Surender: Led the main data sourcing efforts, conducted literature reviews, and assisted with code development and testing.

- Kate Barouch: Led the main coding efforts, contributed to  testing and debugging, reviewed code documentation for JBrowse, and assisted with literature review.
  
- Mariette Peutz: Led literature review, contributed to data sourcing efforts, and assisted with with code development and testing.

---
