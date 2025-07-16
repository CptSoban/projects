
from Bio import SeqIO
import requests
import re
import os
import time

# Define a function to fetch the latest InterProScan version
# def latest_interpro_version():
#     # Fetch the latest version number
#     url = "https://github.com/ebi-pf-team/interproscan"
#     retries = 0
#     max_retries = 5
#     # Retry the request if a ConnectionError occurs
#     while retries < max_retries:
#         try:
#             response = requests.get(url)
#             break
#         except requests.exceptions.ConnectionError as error:
#             print(f"ConnectionError occurred. Retrying {retries}/{max_retries} in 5 seconds...")
#             retries += 1
#             time.sleep(5)

#     # Parse the response to extract the version number    
#     if response.status_code == 200:
#         match = re.search(r'href="/ebi-pf-team/interproscan/releases/tag/([0-9]\.[0-9].\-[0-9]*\.[0-9])', response.text)
#         if match:
#             version = match.group(1)

#             return version

#         else:
#             print("Could not find the version number.")
#     else:
#         print("Failed to fetch the page.")

# Define a function to download the latest InterProScan data
# def get_latest_interpro_data():        
#     # Construct the download URL
#     version = latest_interpro_version()
#     tarball_path = f"{os.getcwd()}/resources/interproscan-data-{version}.tar.gz"
#     data_path = f"{os.getcwd()}/resources/interproscan-{version}/data"
#     # Check if the data directory already exists
#     if os.path.isdir(f"{data_path}"):
#         return data_path
#     else:
#         # Construct the download URL
#         download_url = f"http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{version}/alt/interproscan-data-{version}.tar.gz"
#         print(f"Downloading newest InterProScan version: {download_url}")
        
#         # Download the tarball and extract the data
#         curl_command = f"curl -O --output-dir resources/ {download_url}"
#         print(f"Decompressing InterProScan data: {tarball_path}")
#         tar_command = f"tar -pxzf {tarball_path} -C resources/"
#         remove_tat_command = f"rm -r {tarball_path}"
#         os.system(curl_command)
#         os.system(tar_command)
#         return data_path

# latest_interpro_data = get_latest_interpro_data()

# Define a function to filter sequences that end with an asterisk
def filter_sequences_with_asterisk(input_file, output_file):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        # Parse the input fasta file
        sequences = SeqIO.parse(input_handle, "fasta")

        # Filter sequences that end with an asterisk and do not contain an asterisk elsewhere
        filtered_sequences = (seq for seq in sequences if seq.seq.endswith("*") and seq.seq.count("*") == 1 or "*" not in seq.seq)

        # Remove the trailing asterisk from these sequences
        cleaned_sequences = (seq[:len(seq.seq) - 1] + seq.seq[-1:].replace("*", "") for seq in filtered_sequences)

        # Write the filtered sequences to the output fasta file
        SeqIO.write(cleaned_sequences, output_handle, "fasta")

    return output_file

# Define a rule to filter sequences that end with an asterisk
rule filter_sequences_with_asterisk:
    input:
        all_protein_sequences = config["protein_fasta"]
    output:
        cleaned_aa = "results/{run_id}/functional_annotations/cleaned.aa", 
    run:
        filter_sequences_with_asterisk(input[0], output[0])

# Pull the latest InterProScan Singularity container
rule pull_interpro_sif:
    output:
        f"resources/interproscan_5.75-106.0.sif"
    shell:
        "apptainer pull --force --dir resources/ docker://interpro/interproscan:5.75-106.0"

# Run InterProScan
rule interproscan_run:
    input:
        cleaned_aa = "results/{run_id}/functional_annotations/cleaned.aa",
        interpro_data = "resources/interproscan-5.75-106.0/data",
        container_file = f"resources/interproscan_5.75-106.0.sif"
 
    output:
        interpro_gff = "results/{run_id}/functional_annotations/interproscan/{run_id}.gff3",

    params:
        interpro_output = "results/{run_id}/functional_annotations/interproscan/{run_id}",
        go_terms = "-goterms" if config["GO_terms"] == "yes" else "",

    threads: config["threads"]

    shell:  """ apptainer exec \
            -B {input.interpro_data}:/opt/interproscan/data/ \
            {input.container_file} \
            /opt/interproscan/interproscan.sh \
            -cpu {threads} \
            {params.go_terms} \
            -i {input.cleaned_aa} \
            -b {params.interpro_output} """

# Create a GFF database from the InterProScan results
rule interpro_db:
    input:
        interpro_gff = "results/{run_id}/functional_annotations/interproscan/{run_id}.gff3",

    output:
        interpro_results_db = "resources/{run_id}/interpro_results_db",
    
    conda:
        "../envs/gffutils_db.yaml"
    
    script:
        "../scripts/interproscan_db.py"
