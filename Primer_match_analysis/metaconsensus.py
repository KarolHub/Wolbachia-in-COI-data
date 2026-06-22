#!/usr/bin/env python3

import sys
from collections import Counter, defaultdict
from Bio import SeqIO
import os
from datetime import datetime

# IUPAC codes for ambiguity
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

def read_fasta(file_path):
    sequences = [str(record.seq).upper() for record in SeqIO.parse(file_path, "fasta")]
    return sequences

def validate_sequences(files):
    length = None
    for f in files:
        seqs = read_fasta(f)
        for seq in seqs:
            if length is None:
                length = len(seq)
            elif len(seq) != length:
                raise ValueError(f"Sequences in {f} are not of the same length.")
            if any(base not in "ACGT" for base in seq):
                raise ValueError(f"File {f} contains ambiguous bases. Only A, C, G, T are allowed.")
    return length

def compute_weighted_counts(files):
    file_weights = []
    file_counts = []

    for file in files:
        sequences = read_fasta(file)
        weight = 1.0 / len(sequences)
        file_weights.append(weight)

        position_counts = [Counter() for _ in range(len(sequences[0]))]

        for seq in sequences:
            for i, base in enumerate(seq):
                position_counts[i][base] += weight

        file_counts.append(position_counts)

    return file_counts, len(files)

def compute_metaconsensus(file_counts, num_files, threshold):
    seq_length = len(file_counts[0])
    consensus = []
    frequency_table = {'A': [], 'C': [], 'G': [], 'T': []}

    for i in range(seq_length):
        total_counts = defaultdict(float)
        for file_count in file_counts:
            for base, value in file_count[i].items():
                total_counts[base] += value

        # Normalize by number of files and convert to percentages
        for base in "ACGT":
            freq = (total_counts.get(base, 0.0) / num_files) * 100
            frequency_table[base].append(round(freq, 2))

        passing_bases = {base for base in "ACGT" if frequency_table[base][-1] >= threshold}

        if passing_bases:
            consensus_base = IUPAC_CODES[frozenset(passing_bases)]
        else:
            consensus_base = 'N'

        consensus.append(consensus_base)

    return ''.join(consensus), frequency_table

def write_frequency_table(frequency_table, output_prefix):
    freq_filename = f"{output_prefix}_frequencies.txt"
    with open(freq_filename, 'w') as f:
        for base in "ACGT":
            freqs = ", ".join(str(val) for val in frequency_table[base])
            f.write(f"{base}: {freqs}\n")
    print(f"Base frequency table written to {freq_filename}")

def main():
    if len(sys.argv) < 3:
        print("Usage: metaconsensus.py <threshold_percent> <fasta1> <fasta2> ...")
        sys.exit(1)

    try:
        threshold = float(sys.argv[1])
        assert 0 < threshold <= 100
    except:
        print("Threshold must be a number between 0 and 100.")
        sys.exit(1)

    files = sys.argv[2:]
    try:
        validate_sequences(files)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)

    file_counts, num_files = compute_weighted_counts(files)
    consensus_seq, frequency_table = compute_metaconsensus(file_counts, num_files, threshold)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename = f"metaconsensus_{timestamp}.fasta"
    with open(output_filename, 'w') as f:
        f.write(f">metaconsensus\n{consensus_seq}\n")

    print(f"Metaconsensus written to {output_filename}")

    write_frequency_table(frequency_table, f"metaconsensus_{timestamp}")


if __name__ == "__main__":
    main()
