#!/bin/bash

'''
INPUT FILE: .fastq or .fasta.gz
OUTPUT FILE: .fasta 
'''

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 input_file.fastq[.gz]"
    exit 1
fi

# Input file
input_file="$1"

# Output file (.fasta file)
output_file="${input_file%.fastq*}.fasta"

# Check if input file is gzipped
if [[ "$input_file" == *.gz ]]; then
    # Unzip if it is gzipped
    gunzip -c "$input_file" > "${input_file%.gz}"
    input_file="${input_file%.gz}" # .fastaq file   
fi

# Convert fastq to fasta
sed -n '1~4s/^@/>/p;2~4p' "$input_file" > "$output_file"

# Remove input file (.fastaq file)
rm "$input_file"

echo "Conversion complete. Output saved to $output_file"

