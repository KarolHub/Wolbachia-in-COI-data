#! /usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import IUPACData
import sys

# Extend IUPAC to include Inosine (I) that matches A, C, G, T
IUPACData.ambiguous_dna_values['I'] = {'A', 'C', 'G', 'T'}

# Function to check if two sequences match, considering IUPAC codes and gaps
def match_sequences(primer_seq, reference_seq):
    if len(primer_seq) != len(reference_seq):
        return False  # Sequences must be of the same length

    for p_base, r_base in zip(primer_seq, reference_seq):
        # If there's a gap in either sequence, consider it a match
        if p_base == '-' or r_base == '-':
            continue
        # If the bases are equal, continue checking
        if p_base == r_base:
            continue
        # If the primer base is an IUPAC code, check if it matches the reference base
        if p_base in IUPACData.ambiguous_dna_values:
            if r_base not in IUPACData.ambiguous_dna_values[p_base]:
                return False
        else:
            return False
    return True

# Read in the primers and references from the FASTA files
def process_sequences(primers_file, references_file):
    primers = list(SeqIO.parse(primers_file, "fasta"))
    references = list(SeqIO.parse(references_file, "fasta"))

    total_references = len(references)
    results = []

    for primer in primers:
        primer_name = primer.id
        primer_seq = str(primer.seq)

        matched_references = []
        matched_count = 0

        # Compare the primer to each reference
        for reference in references:
            reference_name = reference.id
            reference_seq = str(reference.seq)

            if match_sequences(primer_seq, reference_seq):
                matched_references.append(reference_name)
                matched_count += 1

        # Calculate percentage match
        if total_references > 0:
            percentage_match = (matched_count / total_references) * 100
        else:
            percentage_match = 0

        # Store the result
        results.append({
            "primer_name": primer_name,
            "percentage_match": round(percentage_match, 2),
            "matched_count": matched_count,
            "matched_references": ','.join(matched_references)
        })

    return results

# Write the results to a file
def write_results(results, output_file):
    with open(output_file, 'w') as f:
        for result in results:
            f.write(f"{result['primer_name']}\t{result['percentage_match']}\t{result['matched_count']}\t{result['matched_references']}\n")

# Main function to run the script
def main():
    primers_file = "primers.fasta"
    references_file = "references.fasta"
    output_file = "output.txt"

    # Process the sequences
    results = process_sequences(primers_file, references_file)

    # Write the results to the output file
    write_results(results, output_file)

    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()