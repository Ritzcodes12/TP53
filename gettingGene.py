from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# Set your email address for NCBI
Entrez.email = "ritojaparia@gmail.com"

# Function to fetch the TP53 mRNA, CDS, and gene sequence
def fetch_tp53_sequences():
    # Fetching the full gene sequence (including introns) from the Gene database
    print("\nFetching full gene sequence...")
    handle = Entrez.efetch(db="nucleotide", id="NC_000017.11", seq_start=7668402, seq_stop=7687550, rettype="gb", retmode="text")
    gene_record = SeqIO.read(handle, "gb")
    handle.close()
    
    print("\nFull gene sequence (includes exons and introns):")
    print(gene_record.seq)
    
    # NCBI identifier for TP53 (Homo sapiens)
    gene_id = "NM_000546"  # This is a RefSeq ID for the TP53 mRNA sequence

    # Fetch the mRNA sequence from NCBI
    print("Fetching mRNA transcript...")
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()
    
    # Print the full DNA sequence (mRNA in DNA form)
    print("mRNA transcript sequence (DNA):")
    print(record.seq)
    
    # Convert the DNA sequence to RNA (replace T with U)
    mrna = record.seq.transcribe()
    print("\nmRNA transcript sequence (RNA):")
    print(mrna)
    
    # Extract the CDS (coding sequence)
    for feature in record.features:
        if feature.type == "CDS":
            cds = feature.extract(record.seq)
            cds_rna = cds.transcribe()
            cds_protein = cds.translate(to_stop=True)
            print("\nCDS sequence (DNA):")
            print(cds)
            print("\nCDS sequence (RNA):")
            print(cds_rna)
            print("\nProtein sequence:")
            print(cds_protein)
            # Write the protein sequence to a FASTA file
            with open("GeneSequences.fasta", "w") as fasta_file:
                SeqIO.write(gene_record, fasta_file, "fasta")
                print("FASTA file written successfully.")

# Call the function
fetch_tp53_sequences()
