import argparse
import pandas as pd

def csv_to_fasta(csv_file, output_fasta_file):

    """ Function to convert a CSV file to a FASTA file """
    # The FASTA file is stored at the same directoey as the csv file
    # The FASTA file has the same name (different extension) as the csv file

    df = pd.read_csv(csv_file) # Read the CSV file into a DataFrame
    with open(output_fasta_file, 'w') as fasta_file:
        for index, row in df.iterrows(): # Iterate through each row in the DataFrame          
            header = f'>{index + 1}' # Create the FASTA header using the row index
            sequence = ''.join([str(int(val)) if val.is_integer() else str(val) for val in row.values]) # Concatenate the values in the row to form the sequence
            fasta_file.write(f'{header}\n{sequence}\n') # Write the header and sequence to the FASTA file


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Path to a CSV file.')
    parser.add_argument('csv_file', type=str, help='The path to the input file (either CSV or FASTA).')
    args = parser.parse_args()

    output_fasta_file_path = args.csv_file.rsplit('.', 1)[0] + '.fa' # Generate the output FASTA file path by replacing the extension
    csv_to_fasta(args.csv_file, output_fasta_file_path) # Convert the CSV file to a FASTA file