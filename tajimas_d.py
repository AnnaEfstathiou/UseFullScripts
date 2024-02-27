import argparse
import pandas as pd
import subprocess
from Bio import SeqIO
import os

def csv_to_fasta(csv_file_path, output_fasta_file_path):

    """ Function to convert a CSV file to a FASTA file """

    df = pd.read_csv(csv_file_path) # Read the CSV file into a DataFrame
    with open(output_fasta_file_path, 'w') as fasta_file:
        for index, row in df.iterrows(): # Iterate through each row in the DataFrame          
            header = f'>{index + 1}' # Create the FASTA header using the row index
            sequence = ''.join([str(int(val)) if val.is_integer() else str(val) for val in row.values]) # Concatenate the values in the row to form the sequence
            fasta_file.write(f'{header}\n{sequence}\n') # Write the header and sequence to the FASTA file


def preprocess_fasta(input_fasta, output_fasta):

    """ Preprocesses the FASTA file to remove lines where all positions are 0. """

    # Open input and output files
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        # Parse the input FASTA file
        for record in SeqIO.parse(input_handle, "fasta"):
            # Check if the sequence contains only zeros
            if '1' in record.seq:
                # Write the sequence to the output file
                SeqIO.write(record, output_handle, "fasta")


def run_tajimas_d(fasta_file_path):

    """ Function to run Tajima's D on a FASTA file """

    safe_fasta_file_path = f"'{fasta_file_path}'" # Ensure the file path is safely quoted to prevent shell interpretation issues
    tajimas_D_command = f'tajimas_d -f {safe_fasta_file_path} -p -t -w' # Construct the command to run Tajima's D
    # subprocess.run(tajimas_D_command, shell=True, check=True) # Execute the command in the shell
    result = subprocess.run(tajimas_D_command, shell=True, capture_output=True, text=True) # Execute the command in the shell and capture output
    # Parsing the results to create a dict
    parsed_result = {}
    lines = result.stdout.split('\n')
    for line in lines:
        if line.strip(): # ignoring empty lines
            key, value = line.split(':')
            parsed_result[key.strip()] = value.strip()
    return parsed_result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a CSV file to a FASTA file or run Tajima\'s D on a FASTA file.')
    parser.add_argument('input_file_path', type=str, help='The path to the input file (either CSV or FASTA).')
    args = parser.parse_args()

    file_path = args.input_file_path # Extract the input file path from the parsed arguments
    file_extension = os.path.splitext(file_path)[1] # Get the file extension of the input file

    if file_extension.lower() == '.csv': # Check if the input file is a CSV file
        output_fasta_file_path = file_path.rsplit('.', 1)[0] + '.fa' # Generate the output FASTA file path by replacing the extension
        csv_to_fasta(file_path, output_fasta_file_path) # Convert the CSV file to a FASTA file
        print(f'Converted CSV to FASTA.')

        preprocessed_file_path = file_path.rsplit('.', 1)[0] + '_processed.fa'
        preprocess_fasta(output_fasta_file_path, preprocessed_file_path)
        print(f'Processed FASTA file (removed all 0s seqs).')

        results = run_tajimas_d(preprocessed_file_path) # Run Tajima's D on the generated FASTA file
        
        rm_command = f'rm "{output_fasta_file_path}" "{preprocessed_file_path}"' # Construct the command to delete the generated FASTA file
        subprocess.run(rm_command, shell=True, check=True) # Execute the command in the shell

    elif file_extension.lower() == '.fasta' or file_extension.lower() == '.fa':

        preprocessed_file_path = file_path.rsplit('.', 1)[0] + '_processed.fa'
        preprocess_fasta(file_path, preprocessed_file_path)
        print(f'Processed FASTA file (removed all 0s seqs).')

        results = run_tajimas_d(preprocessed_file_path)

        rm_command = f'rm "{preprocessed_file_path}"' # Construct the command to delete the generated FASTA file
        subprocess.run(rm_command, shell=True, check=True) # Execute the command in the shell
    
    else:
        raise ValueError('The input file must be either a CSV or a FASTA file.')

    # print("Results:", results) # Print results as a dictionary
