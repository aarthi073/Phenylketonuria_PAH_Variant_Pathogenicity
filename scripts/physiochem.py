physchem = []
with open("../docs/substitutions.txt", "r") as f:
        for line in f:
                clean_line = line.rstrip('\n')
                mut = clean_line[-1]
                if ('G' in mut) | ('A' in mut) | ('V' in mut) | ('L' in mut) | ('I' in mut) | ('P' in mut) | ('M' in clean_line): 
                        hpbal = "hydrophobic_aliphatic"
                        physchem.append(hpbal)
                elif ('F' in mut) | ('Y' in mut) | ('W' in mut):
                        hpbar = "hydrophobic_aromatic"
                        physchem.append(hpbar)
                elif ('S' in mut) | ('T' in mut) | ('C' in mut) | ('N' in mut) | ('Q' in mut):
                        hphilP = "hydrophilic_polar"
                        physchem.append(hphilP)                                      
                elif ('D' in mut) | ('E' in mut):
                        hphilA = "hydrophilic_acidic"
                        physchem.append(hphilA)
                elif ('R' in mut) | ('H' in mut) | ('K' in mut):
                        hphilB = "hydrophilic_basic"
                        physchem.append(hphilB)
with open(snakemake.output.txt_file, "w") as f:
        for i in physchem:
                f.write(i+"\n")
