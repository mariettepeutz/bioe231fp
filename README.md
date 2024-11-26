# The Flavi Fighters
## BioE231 Final Project by Smrithi Surender, Marriëtte Peutz, & Kate Barouch

This project provides a customizable database installer for JBrowse2, focused on the Flavivirus family, to support pandemic researchers by integrating genome alignments, annotations, and bioinformatics tools for rapid pathogen identification and evolutionary analysis.

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

Additional Notes:
- If any errors occur, the script will exit and display an error message.

## Flavivirus Family

We chose to focus on the Flavivirus family due to the increasing global health threats posed by this virus family with climate warming expanding mosquito habitats and contributing to outbreaks in new regions.

### Dengue Virus (DENV): 
A major global health threat, with over 6,800 cases reported in the United States as of October 31, 2024—more than double the cases reported in 2023—driven by climate change and increasing prevalence in the Americas and Asia.

### Zika Virus (ZIKV): 
Recognized for causing severe congenital outcomes like microcephaly, it is included on the WHO priority list for diseases with epidemic potential, and there is currently no vaccine available.

### West Nile Virus (WNV): 
A widespread virus transmitted by Culex mosquitoes, with severe neurological complications in some cases, and no available vaccine, underscoring the importance of further research.

# Need to finish writing this section as we go!
## Table of Contents
### x
### x

## Installation Instructions:
### Prerequisites - Any required software or libraries (e.g., JBrowse2, AWS setup).
### Database Installer Setup - Step-by-step instructions on how to download and set up the database installer.
### Running the Database - Instructions to launch the JBrowse2 instance and customize configurations.

## Data Sources - Describe data sources used, with clear distinctions between external sources and any self-generated data.

## Usage Guide:
### Navigating JBrowse2 - Explain navigation options, annotation tracks, and visualization tools.
### Customizing Data Views - Instructions for users to customize the data (e.g., adding new viral strains).

## Collaborator Contributions - Brief statements on each team member’s contributions.

