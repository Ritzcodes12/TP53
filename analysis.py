import json

# Load the JSON data from a file
with open("mutation_data.json", "r") as json_file:
    data = json.load(json_file)

# Open a file to write the mutation data
with open("mutation_data.txt", "w") as file:
    for mutation in data:
        # Extract relevant fields
        uuid = mutation['uuid']
        start_position = mutation['mutation_info']['data']['start_position']
        end_position = mutation['mutation_info']['data']['end_position']
        gene_aa_change = ", ".join(mutation['mutation_info']['data']['gene_aa_change'])
        genomic_dna_change = mutation['mutation_info']['data']['genomic_dna_change']
        cosmic_ids = ", ".join(mutation['mutation_info']['data'].get('cosmic_id', []))

        # Write the extracted information to the file
        file.write(f"UUID: {uuid}\n")
        file.write(f"Start Position: {start_position}\n")
        file.write(f"End Position: {end_position}\n")
        file.write(f"Gene AA Change: {gene_aa_change}\n")
        file.write(f"Genomic DNA Change: {genomic_dna_change}\n")
        file.write(f"COSMIC IDs: {cosmic_ids}\n")  # Include COSMIC IDs if available
        file.write("\n")  # Add a blank line for better readability

print("Mutation data has been saved to mutation_data.txt.")
