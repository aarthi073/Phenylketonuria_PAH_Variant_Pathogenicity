from Bio import Entrez
import json
import requests

#Entrez is a global search system by NCBI to search multiple health and science databases 

search = Entrez.esearch(db="clinvar", term="PAH AND [p.Asp229Gly, p.Gly239Ala, p.Phe263Ser, p.Ala342Pro]")

record = Entrez.read(search)


#ids = record["idlist"]

#print(record)


file = "./esearch.json"

with open(file, "r") as f:

        file_dict = json.load(f)
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

#NCBI stores results on server and returns tokens for retrieval; WebEnv behaves like a session ID for each E-utility call, and query_key is an integer identification of which result within the session is needed.
        idlist = file_dict["esearchresult"]["idlist"]

        id_string = ",".join(idlist)


        params = {
	    "db": "clinvar",
	    "id": id_string,
	    "retmode": "json"
	}

print(idlist)


response = requests.post(url, data=params)
response.raise_for_status()
ids = response.json()

with open("esummary.json", "w") as new_f:
        json.dump(ids, new_f, indent=2)


