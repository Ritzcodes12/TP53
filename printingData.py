import pandas as pd
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Load the reference TP53 sequence (FASTA file)
reference_seq = SeqIO.read('Gene_sequences.fasta', 'fasta')
ref_seq = reference_seq.seq

print(f"Reference sequence length: {len(ref_seq)}")

# Define the starting position of the reference sequence
reference_start_position = 7668401  # Adjust as per your reference genome

# Load the TSV data into a pandas DataFrame
mutations_df = pd.read_csv('TP53_mutations.tsv', sep='\t')

def apply_mutation(ref_seq, dna_change, protein_change):
    print(f"{i}: Processing mutation: {dna_change}")  
    chrom, change = dna_change.split(':g.')
    mutated_seq = MutableSeq(str(ref_seq))

    if '>' in change:  # Substitution
        pos, change_type = change.split('>')
        mutation_position = int(pos[:-1])
        ref_nuc = pos[-1]
        alt_nuc = change_type

        mutation_index = mutation_position - reference_start_position - 1

        if 0 <= mutation_index < len(ref_seq):
            if ref_seq[mutation_index] == ref_nuc:
                mutated_seq[mutation_index] = alt_nuc
            else:
                print(f"Reference mismatch at {mutation_position}: expected {ref_nuc}, found {ref_seq[mutation_index]}.")
        else:
            print(f"Mutation index {mutation_position} out of bounds.")

    elif 'del' in change:  #Eg: Deletion chr17:g.7675185delCAGGGCAGGTCTTG
        genomic_position = int(change.split('del')[0])  # Extract the genomic position
        deletion_seq = change.split('del')[1]  # Extract the sequence being deleted

        # Calculate the start index relative to the reference sequence
        start_index = genomic_position - reference_start_position - 1  # Convert to 0-based index
        end_index = start_index + len(deletion_seq) - 1  # Calculate the end index based on deletion length
    
        # Deleting the sequence from the mutated sequence
        if 0 <= start_index <= end_index < len(ref_seq):
            del mutated_seq[start_index:end_index + 1]
        else:
            print(f"Deletion range {start_index}-{end_index} out of bounds.")

    elif 'ins' in change:  # Insertion
        pos, inserted_seq = change.split('ins')
        insertion_position = int(pos.split('_')[0])
        start_index = insertion_position - reference_start_position - 1

        if 0 <= start_index < len(ref_seq):
            mutated_seq = mutated_seq[:start_index + 1] + inserted_seq + mutated_seq[start_index + 1:]
        else:
            print(f"Insertion index {start_index} out of bounds.")
    else:
        print(f"Unknown mutation type in change: {change}")
        return None
    print(mutated_seq)
    return mutated_seq

# Process mutations and save sequences
mutated_sequences = []
i=1
for index, row in mutations_df.iterrows():
    dna_change = row['dna_change']
    protein_change = row['protein_change']
    type1 = row['type']
    try:
        mutated_seq = apply_mutation(ref_seq, dna_change, protein_change)
        i=i+1
        if mutated_seq:
            seq_record = SeqRecord(mutated_seq, id=row['ssm_id'], description=f"{protein_change}, Type: {row['type']}")
            mutated_sequences.append(seq_record)
            print(f"Mutation applied for {row['ssm_id']}: {protein_change} : {dna_change} :{type1}")
        else:
            print(f"Mutation {row['ssm_id']} could not be applied.")
    except (IndexError, ValueError) as e:
        print(f"Error processing mutation {row['ssm_id']}: {e}")

# Save the mutated sequences
if mutated_sequences:
    SeqIO.write(mutated_sequences, 'SEQUENCES.fasta', 'fasta')
    print("Mutated sequences have been saved to 'SEQUENCES.fasta'.")
else:
    print("No mutated sequences to save.")
