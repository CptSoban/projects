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

echo "Removing FASTQ files with fewer than 10 reads..."

# Function: return success if file has at least 10 reads
has_min10_reads() {
    zcat -f -- "$1" 2>/dev/null | awk '
        NR >= 40 { exit 0 }     # 10 reads = 40 lines
        END {
            if (NR >= 40) exit 0
            else exit 1
        }
    '
}

# Loop over renamed files and remove if fewer than 10 reads
for fq in "$output_dir"/*.fastq.gz; do
    [ -e "$fq" ] || continue

    if ! has_min10_reads "$fq"; then
        echo "Removing low-read file (<10 reads): $fq"
        rm "$fq"
    fi
done

echo "Demultiplexing complete."
done

