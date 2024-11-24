from Bio import SeqIO
from Bio.Align import PairwiseAligner

def read_sequences(original_file, mutated_file):
    # Read the original sequence from the original FASTA file
    original_seq = list(SeqIO.parse(original_file, "fasta"))[0].seq
    
    # Read the mutated sequences from the mutated FASTA file
    mutated_seqs = [seq.seq for seq in SeqIO.parse(mutated_file, "fasta")]
    
    return original_seq, mutated_seqs

def align_and_show(original_seq, mutated_seqs):
    # Initialize the pairwise aligner
    aligner = PairwiseAligner()
    
    # Set alignment mode to global (or local if needed)
    aligner.mode = 'global'  # You can change to 'local' for local alignment

    aligner.open_gap_score = -10  # Example penalty for opening a gap
    aligner.extend_gap_score = -0.5  # Example penalty for extending a gap

    
    # Iterate over each mutated sequence and align it with the original sequence
    for i, mutated_seq in enumerate(mutated_seqs):
        alignments = aligner.align(original_seq, mutated_seq)
        
        # Get the best alignment
        best_alignment = alignments[0]  # The first alignment in the list (best one)
        
        # Print the alignment in a human-readable format
        print(f"\nAlignment for Mutated Sequence {i + 1}:")
        print(best_alignment)  # This shows the target, query, and alignment representation
        
        # Print the alignment score
        print(f"Alignment score for Mutated Sequence {i + 1}: {best_alignment.score}")

# File paths to the original and mutated FASTA files
original_file = "Gene_sequences.fasta"
mutated_file = "SEQUENCES.fasta"

# Read the sequences
original_seq, mutated_seqs = read_sequences(original_file, mutated_file)

# Align the sequences and print the alignment and scores
align_and_show(original_seq, mutated_seqs)
