import pandas as pd
import re  # Importing regex for flexible position extraction

# Function to read domain specifications
def load_domains(domain_file):
    # Assuming the domain file has columns like: "Domain/Region", "Start", "End"
    domains_df = pd.read_csv(domain_file, sep="\t")
    domains = []
    for _, row in domains_df.iterrows():
        domains.append({
            "name": row["Domain/Region"],
            "start": float(row["Start"]),
            "end": float(row["End"])
        })
    return domains

# Function to read mutation information and extract positions
def load_mutations(mutation_file):
    # Assuming the mutation file has a column "protein_change"
    mutations_df = pd.read_csv(mutation_file, sep="\t")
    mutations = []
    for _, row in mutations_df.iterrows():
        protein_change = row["protein_change"]
        position = None  # Initialize position as None
        
        try:
            # Using regex to extract numeric part indicating the position
            if re.search(r'R(\d+)', protein_change):  # For mutations like R175H
                position = int(re.search(r'R(\d+)', protein_change).group(1))
            elif re.search(r'V(\d+)', protein_change):  # For mutations like V143Afs*14 or V147_D148du
                position = int(re.search(r'V(\d+)', protein_change).group(1))
            elif re.search(r'X(\d+)', protein_change):  # For splice mutations
                position = int(re.search(r'X(\d+)', protein_change).group(1))
            elif re.search(r'M(\d+)', protein_change):  # For mutations like M66Pfs*60
                position = int(re.search(r'M(\d+)', protein_change).group(1))
            
            # Append only if position is found
            if position is not None:
                mutations.append({
                    "mutation": protein_change,
                    "position": position
                })
        except ValueError:
            # Skip invalid entries
            continue
            
    return mutations

# Function to find affected domains
def find_affected_domains(mutations, domains):
    for mutation in mutations:
        affected = []
        for domain in domains:
            if domain["start"] <= mutation["position"] <= domain["end"]:
                affected.append(domain["name"])
        
        if affected:
            print(f"Mutation: {mutation['mutation']} at Position: {mutation['position']} affects Domain(s): {', '.join(affected)}")
        else:
            print(f"Mutation: {mutation['mutation']} at Position: {mutation['position']} does not affect any known domain")

# File paths (update these paths with your actual file paths)
domain_file = "domain.tsv"
mutation_file = "TP53_mutations.tsv"

# Load domains and mutations
domains = load_domains(domain_file)
mutations = load_mutations(mutation_file)

# Find and print affected domains
find_affected_domains(mutations, domains)
