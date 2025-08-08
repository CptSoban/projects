 #!/bin/bash
 
 # Demultiplexing reads using Nanopore Dorado
 dorado demux ${snakemake_input} \
     -o ${snakemake_output[output_dir]} \
     --emit-fastq \
     --kit-name ${snakemake_params[kit_name]} \
     --threads ${snakemake_threads}

echo "Compressing demultiplexed reads..."

# Compress the demultiplexed reads to save space
gzip ${snakemake_output[output_dir]}/*.fastq

echo "Renaming demultiplexed reads..."

# Rename the demultiplexed reads based on the barcode
output_dir="${snakemake_output[output_dir]}"
for file in "$output_dir"/*.fastq.gz; do
    base_name=$(basename "$file")
    if [[ "$base_name" =~ _barcode[0-9]+\.fastq\.gz$ ]]; then
        new_name=$(echo "$base_name" | sed -E 's/^.*_(barcode[0-9]+\.fastq\.gz)$/\1/')
    elif [[ "$base_name" =~ _unclassified\.fastq\.gz$ ]]; then
        new_name="unclassified.fastq.gz"
    else
        continue
    fi
    
    mv "$file" "$output_dir/$new_name"

echo "Removing empty FASTQ files..."

# Function to check if a FASTQ is empty (0 reads)
is_fastq_empty() {
    # Looks for at least one non-empty sequence line (2nd line of each 4-line FASTQ block)
    zcat "$1" | awk 'NR%4==2 && length($0) > 0 { found=1; exit } END { exit !found }'
}

# Loop over renamed files and remove if empty
for fq in "$output_dir"/*.fastq.gz; do
    if ! is_fastq_empty "$fq"; then
        echo "Removing empty file: $fq"
        rm "$fq"
    fi
done

echo "Demultiplexing complete."
done

