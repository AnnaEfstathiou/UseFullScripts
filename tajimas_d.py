import argparse
import pandas as pd
import subprocess
import os

def csv_to_fasta(csv_file_path, output_fasta_file_path):
    df = pd.read_csv(csv_file_path)
    with open(output_fasta_file_path, 'w') as fasta_file:
        for index, row in df.iterrows():
            header = f'>{index + 1}'
            sequence = ''.join([str(int(val)) if val.is_integer() else str(val) for val in row.values])
            fasta_file.write(f'{header}\n{sequence}\n')

def run_tajimas_d(fasta_file_path):
        # Ensure the file path is safely quoted to prevent shell interpretation issues
    safe_fasta_file_path = f"'{fasta_file_path}'"
    command = f'tajimas_d -f {safe_fasta_file_path} -p -t -w'
    subprocess.run(command, shell=True, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a CSV file to a FASTA file or run Tajima\'s D on a FASTA file.')
    parser.add_argument('input_file_path', type=str, help='The path to the input file (either CSV or FASTA).')
    args = parser.parse_args()

    file_path = args.input_file_path
    file_extension = os.path.splitext(file_path)[1]

    if file_extension.lower() == '.csv':
        output_fasta_file_path = file_path.rsplit('.', 1)[0] + '.fasta'
        csv_to_fasta(file_path, output_fasta_file_path)
        print(f'Converted CSV to FASTA: {output_fasta_file_path}')
        run_tajimas_d(output_fasta_file_path)
    
    elif file_extension.lower() == '.fasta':
        run_tajimas_d(file_path)
    
    else:
        raise ValueError('The input file must be either a CSV or a FASTA file.')

