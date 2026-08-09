import json
import re

file="../docs/esummary.json"


with open(file, "r") as f:
        file_dict = json.load(f)

ids  = file_dict['result']['uids']


vlist = []
mut_sub=[]
mut_classification = []
var_type = []
#Extract the variant name, HGVS protein notation, clinical significance label
for uid in ids:
        record = file_dict['result'][uid]
        #use .get() to prevent crashes from brackets if there is an empty field
        variant = record.get('obj_type', 'Unknown')
        vname = record.get('variation_name', 'Unknown')
        HGVS = record.get('title', 'Unknown')
        mutation_type = record.get('molecular_consequence_list', 'Unknown')
        germ = record.get('germline_classification', 'Unknown')
        clinical = record.get('clinical_impact_classification', 'Unknown')
        onco = record.get('oncogenicity_classification', 'Unknown')
        vlist.append([variant, vname, HGVS, mutation_type, germ, clinical, onco])

#only include missense vairants in vlist
for i in vlist:
        vlist = [i for i in vlist if "missense variant" in i[4]]
for i in vlist:
        desc = i[5]["description"]
        sub=i[2].rpartition(".")[-1]
        if ("?" not in sub) & ("_" not in sub) & ("-" not in sub) & ("+" not in sub) & (">" not in sub):
                mut = sub
        clean_mis = mut.rpartition(")")[0]
        mut_type = i[4]
        for i in mut_type:
                var_type.append(i)
        mut_classification.append(desc)
        mut_sub.append(clean_mis)

from Bio.Data.IUPACData import protein_letters_3to1

mut_1_code = []
for mutation in mut_sub:
        if not mutation.strip():
                continue
        #find all conditions where there are exactly two lowercase letters following an uppercase letter
        original = mutation[:3]
        new = mutation[-3:]

        og_code = protein_letters_3to1[original]
        new_code = protein_letters_3to1[new]

        new_mut = f"{og_code}{mutation[3:-3]}{new_code}"

        mut_1_code.append(new_mut)



join = " ".join(mut_1_code)
print(join)

with open("../docs/cleaned.json","w") as f:
         json.dump(vlist, f, indent=2)

with open("../docs/substitutions.txt", "w") as f:
         for i in mut_1_code:
                f.write(i+"\n")
with open("../docs/label.txt", "w") as f:
         for i in mut_classification:
                f.write(i+"\n")
with open("../docs/variant_type", "w") as f:
         for i in var_type:
                f.write(i+"\n")

