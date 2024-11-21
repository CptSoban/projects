from Bio import SeqIO
import requests
import re
import os
import time

def latest_interpro_version():
    # Fetch the latest version number
    url = "https://github.com/ebi-pf-team/interproscan"
    retries = 0
    max_retries = 5
    while retries < max_retries:
        try:
            response = requests.get(url)
            break
        except requests.exceptions.ConnectionError as error:
            print(f"ConnectionError occurred. Retrying {retries}/{max_retries} in 5 seconds...")
            retries += 1
            time.sleep(5)
            
    if response.status_code == 200:
        match = re.search(r'href="/ebi-pf-team/interproscan/releases/tag/([0-9]\.[0-9].\-[0-9]*\.[0-9])', response.text)
        if match:
            version = match.group(1)

            return version

        else:
            print("Could not find the version number.")
    else:
        print("Failed to fetch the page.")


def get_latest_interpro_data():        
    # Construct the download URL
    version = latest_interpro_version()
    tarball_path = f"{os.getcwd()}/resources/interproscan-data-{version}.tar.gz"
    data_path = f"{os.getcwd()}/resources/interproscan-{version}/data"
    if os.path.isdir(f"{data_path}"):
        return data_path
    else:
        download_url = f"http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/{version}/alt/interproscan-data-{version}.tar.gz"
        print(f"Download URL: {download_url}")
        
        # Execute the curl command to download the file  
        curl_command = f"curl -O --output-dir resources/ {download_url}"
        tar_command = f"tar -pxzf {tarball_path} -C resources/"
        pull_command = f"apptainer pull --dir rna_seq_workflow/resources docker://interpro/interproscan:latest"
        os.system(curl_command)
        os.system(tar_command)
        os.system(pull_command)
        return data_path

latest_interpro_data = get_latest_interpro_data()

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

#If the reference genome was annotated using BRAKER, then there will be asterisks at the of the sequences in .codingseq and .aa
#indicating the end of the specific sequences. These can be removed without concern, which is necessary for Interproscan.
rule filter_sequences_with_asterisk:
    input:
        all_protein_sequences = config["protein_fasta"] #"results/DEG_analysis/{run_id}/{contrast}/{contrast}.aa",
    output:
        cleaned_aa = "results/functional_annotations/{run_id}/cleaned.aa", #"results/functional_annotations/{run_id}/{contrast}/{contrast}_clean.aa",
    run:
        filter_sequences_with_asterisk(input[0], output[0])
        

rule pull_interpro_sif:
    output:
        "resources/interproscan_latest.sif"
    shell:
        "apptainer pull --dir resources/ interproscan_latest.sif docker://interpro/interproscan:latest"


rule interproscan_run:
    input:
        cleaned_aa = "results/functional_annotations/{run_id}/cleaned.aa", #"results/functional_annotations/{run_id}/{contrast}/{contrast}_clean.aa",
        interpro_data = latest_interpro_data

    output:
        interpro_gff = "results/functional_annotations/{run_id}/interproscan/{contrast}.gff3",

    params:
        interpro_output = "results/functional_annotations/{run_id}/interproscan/{contrast}",
        go_terms = "-goterms" if config["GO_terms"] == "yes" else "",

    threads:
        config["threads"]

    shell:  """ apptainer exec \
            -B {input.interpro_data}:/opt/interproscan/data/ \
            resources/interproscan_latest.sif \
            /opt/interproscan/interproscan.sh \
            -cpu {threads} \
            {params.go_terms} \
            -i {input.cleaned_aa} \
            -b {params.interpro_output} """

rule interpro_db:
    input:
        interpro_gff = "results/functional_annotations/{run_id}/interproscan/{contrast}.gff3",

    output:
        interpro_results_db = "resources/{run_id}/{contrast}_interpro_results_db",
    
    conda:
        "../envs/gffutils_db.yaml"
    
    script:
        "../scripts/interproscan_db.py"
