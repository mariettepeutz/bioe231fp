import subprocess
from BCBio import GFF
from Bio import SeqIO
from pathlib import Path

# Our big questions for Bear: 
# How do we get around manualing putting in each URL? 
# Should this somehow connect to the terminal? I'm confused how to create the desired deliverable. 

# List of GenBank URLs to process 
genbank_urls = [
    "https://www.ncbi.nlm.nih.gov/nuccore/OM513935.1?report=genbank&format=text",
    # Add more URLs here
]

# Directory to store downloaded and processed files
output_dir = Path("genbank_files")
output_dir.mkdir(exist_ok=True)

for url in genbank_urls:
    # Extract name each file from the URL
    identifier = url.split("/")[-1].split("?")[0] # this gets us the unique GenBank identifier (ie https://www.ncbi.nlm.nih.gov/nuccore/OM513935.1?report=genbank)
    genbank_filename = output_dir / f"{identifier}.gb"
    gff_filename = output_dir / f"{identifier}.gff3"
    compressed_gff_filename = gff_filename.with_suffix(".gff3.gz")

    # Download the GenBank file
    print(f"Downloading {identifier}...")
    subprocess.run(["wget", url, "-O", str(genbank_filename)])

    # Convert to GFF3
    print(f"Converting {genbank_filename} to GFF3 format...")
    with genbank_filename.open("r") as input_handle, gff_filename.open("w") as output_handle:
        GFF.write(SeqIO.parse(input_handle, "genbank"), output_handle)

    # Compress the GFF3 file with bgzip
    print(f"Compressing {gff_filename}...")
    subprocess.run(["bgzip", str(gff_filename)])

    # Index the compressed GFF3 file with tabix
    print(f"Indexing {compressed_gff_filename}...")
    subprocess.run(["tabix", "-p", "gff", str(compressed_gff_filename)])

    print(f"Completed processing for {identifier}\n")

print("All files processed.")

    ### Next step load annotation track into jbrowse ###






