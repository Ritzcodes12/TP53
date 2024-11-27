import requests
import json

# Define the correct endpoint for the GDC API
url = "https://api.gdc.cancer.gov/ssms"

# Define the JSON filter to find mutations for TP53 (using its Ensembl gene ID)
filters = {
    "op": "in",
    "content": {
        "field": "genes.gene_id",
        "value": ["ENSG00000141510"]  # TP53 Ensembl gene ID
    }
}

# Define the parameters for the request, including size to specify the number of results
params = {
    "filters": json.dumps(filters),  # Convert the filter to a JSON string
    "size": 1000,  # Adjust the size based on how many results you want
    "fields": "ssm_id, gene.gene_id, mutation_id, consequence_type",  # Retrieve mutation-specific fields
    "pretty": "true"
}

# Make the POST request to the GDC API
response = requests.post(url, headers={'Content-Type': 'application/json'}, data=json.dumps(params))

# Check if the request was successful
if response.status_code == 200:
    # Parse the JSON response
    data = response.json()
    # Extract the list of ssm_ids
    mutations = data["data"]["hits"]
    
    print(f"Retrieved {len(mutations)} mutations for TP53:")
    
    # Open a file to write UUIDs
    with open("uuid_list.txt", "w") as file:
        for mutation in mutations:
            uuid = mutation['ssm_id']
            print(uuid)  # Print the UUID to the console
            file.write(f"{uuid}\n")  # Write UUID to the file
else:
    print(f"Error: {response.status_code}")
    print(response.text)
