#!/usr/bin/env python3

import sys
from collections import Counter
from Bio import SeqIO
import os

# IUPAC nucleotide codes for ambiguity
IUPAC_CODES = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

def get_consensus(sequences, threshold):
    seq_length = len(sequences[0])
    consensus = []

    for i in range(seq_length):
        column = [seq[i] for seq in sequences]
        counts = Counter(column)
        total = sum(counts[nuc] for nuc in "ACGT")
        included = set()

        for nuc in "ACGT":
            if (counts[nuc] / total) * 100 >= threshold:
                included.add(nuc)

        if included:
            consensus.append(IUPAC_CODES[frozenset(included)])
        else:
            consensus.append('N')  # Default if nothing meets threshold

    return ''.join(consensus)

def main():
    if len(sys.argv) != 3:
        print("Usage: consensus.py <threshold_percent> <input_fasta>")
        sys.exit(1)

    threshold = float(sys.argv[1])
    input_file = sys.argv[2]

    sequences = [str(record.seq).upper() for record in SeqIO.parse(input_file, "fasta")]
    if not sequences:
        print("No sequences found in the input FASTA file.")
        sys.exit(1)

    if len(set(len(seq) for seq in sequences)) != 1:
        print("All sequences must be of the same length.")
        sys.exit(1)

    consensus_seq = get_consensus(sequences, threshold)

    output_filename = f"consensus_{os.path.basename(input_file)}"
    with open(output_filename, 'w') as out:
        out.write(f">consensus\n{consensus_seq}\n")

    print(f"Consensus written to {output_filename}")

if __name__ == "__main__":
    main()
