import requests
import json


url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
params = {
    "db": "clinvar",
    "term": 'PAH[gene] AND Phenylketonuria[disease/phenotype] AND (pathogenic[Clinical significance] OR likely pathogenic[Clinical significance])',
    "retmax": 500,
    "retmode": "json",
}
response = requests.get(url, params=params)
data = response.json()

with open(snakemake.output.json_file, "w") as f:
    json.dump(data, f, indent=2)

print(data["esearchresult"]["idlist"])

print(data)
