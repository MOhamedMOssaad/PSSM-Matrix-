import random
import math

def calculate_hamming_distance(seq1, seq2):
    return sum(n1 != n2 for n1, n2 in zip(seq1, seq2))

def calculate_total_distance(pattern, sequences, pattern_length):
    total_distance = 0
    for sequence in sequences:
        min_hamming_distance = math.inf
        for i in range(len(sequence) - pattern_length + 1):
            current_distance = calculate_hamming_distance(pattern, sequence[i:i+pattern_length])
            min_hamming_distance = min(min_hamming_distance, current_distance)
        total_distance += min_hamming_distance
    return total_distance

def generate_patterns(pattern_length):
    nucleotides = ['A', 'C', 'G', 'T']
    patterns = []
    generate_patterns_helper(nucleotides, pattern_length, "", patterns)
    return patterns

def generate_patterns_helper(nucleotides, pattern_length, current_pattern, patterns):
    if len(current_pattern) == pattern_length:
        patterns.append(current_pattern)
    else:
        for nucleotide in nucleotides:
            generate_patterns_helper(nucleotides, pattern_length, current_pattern + nucleotide, patterns)

def find_consensus_pattern(dna_sequences, pattern_length):
    sequence_count = len(dna_sequences)
    best_distance = math.inf
    best_pattern = ""

    patterns = generate_patterns(pattern_length)

    for pattern in patterns:
        distance = calculate_total_distance(pattern, dna_sequences, pattern_length)

        if distance < best_distance:
            best_distance = distance
            best_pattern = pattern
    alignment_indices = []
    for index, sequence in enumerate(dna_sequences):
        min_hamming_distance = math.inf
        min_index = 0

        for i in range(len(sequence) - pattern_length + 1):
            current_distance = calculate_hamming_distance(best_pattern, sequence[i:i + pattern_length])
            if current_distance < min_hamming_distance:
                min_hamming_distance = current_distance
                min_index = i

        alignment_indices.append(min_index)
    for index, value in enumerate(alignment_indices):
        print(f"Sequence {index + 1}>>>> {value}")

    return best_pattern

def read_dna_from_file():
    file_path = "rawDNA.txt"
    dna = []
    try:
        with open(file_path, 'r') as f:
            sequence_count, sequence_length, pattern_length = map(int, f.readline().split())
            dna = [line.strip().upper() for line in f.readlines()]
    except FileNotFoundError:
        print("File not found. Please make sure 'rawDNA.txt' exists.")
    return sequence_count, sequence_length, pattern_length, dna

def generate_random_dna(sequence_count, sequence_length):
    # Generate random DNA sequences
    nucleotides = ['A', 'C', 'G', 'T']
    sequences = []
    for _ in range(sequence_count):
        sequence = ''.join(random.choices(nucleotides, k=sequence_length))
        sequences.append(sequence)
    return sequences

def perform_multiple_alignment(dna_sequences, pattern_length):
    consensus = find_consensus_pattern(dna_sequences, pattern_length)
    for index, sequence in enumerate(dna_sequences):
        min_distance = float('inf')
        min_index = 0
        for i in range(len(sequence) - pattern_length + 1):
            distance = calculate_hamming_distance(consensus, sequence[i:i + pattern_length])
            if distance < min_distance:
                min_distance = distance
                min_index = i
    return consensus

random.seed(42)

if __name__ == '__main__':
    while True:
        choice = input('Write 1 to read the prepared raw DNA from a file\nWrite 2 for random DNA generator\nWrite 0 to exit\n')
        if choice == '1':
            sequence_count, sequence_length, pattern_length, dna = read_dna_from_file()
            if dna:
                consensus = perform_multiple_alignment(dna, pattern_length)
                print("Consensus pattern:", consensus)
        elif choice == '2':
            sequence_count = int(input('Please enter the number of sequences: '))
            sequence_length = int(input('Please enter the length of each sequence: '))
            pattern_length = int(input('Please enter the length for the pattern to be found: '))
            dna_sequences = generate_random_dna(sequence_count, sequence_length)
            consensus = perform_multiple_alignment(dna_sequences, pattern_length)
            print("Consensus pattern:", consensus)
        elif choice == '0':
            break
        else:
            print('Wrong choice. Please try again.')
