import requests
import json

# Function to fetch mutation data for a given UUID
def fetch_mutation_data(uuid):
    # Base URL for the GDC mutation endpoint
    base_url = f"https://api.gdc.cancer.gov/ssms/{uuid}"
    
    response = requests.get(base_url)

    if response.status_code == 200:
        return response.json()  # Convert the response to JSON
    else:
        print(f"Error fetching data for UUID {uuid}: {response.status_code} - {response.text}")
        return None

# Load UUIDs from the 'uuid_list.txt' file
with open('uuid_list.txt', 'r') as f:
    uuids = [line.strip() for line in f.readlines()]

# Create a list to store mutation data
mutation_data = []

# Loop through each UUID and fetch mutation data
for uuid in uuids:
    mutation_info = fetch_mutation_data(uuid)
    if mutation_info:
        mutation_data.append({
            'uuid': uuid,
            'mutation_info': mutation_info
        })

# Save the mutation data to a JSON file
with open('mutation_data.json', 'w') as outfile:
    json.dump(mutation_data, outfile, indent=2)

print(f"Retrieved data for {len(mutation_data)} UUIDs. Check 'mutation_data.json' for results.")
