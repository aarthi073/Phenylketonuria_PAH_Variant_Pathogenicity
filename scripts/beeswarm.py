import seaborn as sns
import matplotlib.pyplot as plt
import csv
import sys
import os
import numpy as np
import pandas as pd
from collections import defaultdict


#higher PhyloP100 means higher conservation of the gene
def main():
    try: 
         Phd = snakemake.input.phd_snp
         mutpred2 = snakemake.input.mutpred2
          
#        Phd = sys.argv[1]
 #       mutpred2 = sys.argv[2]
    except IndexError:
            print("no pathogenicity prediction table passed")


    def Phd_SNPg_dataframe(Phd_SNPg):
        headers = ["CHROM", "REF", "ALT", "MUT", "CODING", "PREDICTION", "SCORE", "PhyloP100"]
        cols = ["CHROM", "POS", "ID", "REF", "ALT", "MUT", "CODING", "PREDICTION", "SCORE", "FDR", "PhyloP100", "AvgPhyloP100", "transcript", "gene", "strand", "coordinates(gDNA/cDNA/protein)", "region", "info"]
        rows = []
        with open(Phd_SNPg) as f:
            for line in f:
                  if line.startswith("#"):
                        continue                  
                  fields = line.strip().split(None, 17)
                  rows.append(fields)
        df = pd.DataFrame(rows, columns = cols)
        numeric_cols = ["CHROM", "SCORE", "PhyloP100"]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.drop(columns=["POS", "ID", "FDR", "AvgPhyloP100",  "transcript", "gene", "strand", "coordinates(gDNA/cDNA/protein)", "region", "info"])
        df.to_csv("./docs/Phd-Snp_Results.tsv", sep="\t", index=False)
        return df

# I HAVE TO MAKE A TABLE
    def Mutpred2_dataframe(mutpred2):
#        headers = ["ID", "Substitution", "Score", "Molecular Mechanisms", "Motif Information"]
        headers = ["ID", "Substitution", "MutPred2 score", "Molecular mechanisms with Pr >= 0.01 and P < 0.05", "Motif information"]
        cols = ["ID", "Substitution", "MutPred2 score", "Molecular mechanisms with Pr >= 0.01 and P < 0.05", "Motif information", "Remarks"]        
        rows = []
        with open(mutpred2, "r") as f:
            for line in f:
                cleaned_line = line.strip()
                if (
                    not cleaned_line
                    or cleaned_line.startswith(">")
                    or cleaned_line.startswith("#")
                   ):
                     continue    
                if "," in cleaned_line:
                        fields = [item.strip() for item in cleaned_line.split(",", 5)]
                else:
                        fields = cleaned_line.split(None,5)
                
                rows.append(fields)
        df = pd.DataFrame(rows, columns = cols)
        df = df.drop(columns=["Remarks"])
        df.columns = df.columns.str.strip()        

        #Convert score to numeric to evaluate pathogenicity
        df["MutPred2 score"] = pd.to_numeric(df["MutPred2 score"], errors="coerce")
        df["Pathogenicity"] = np.where(df["MutPred2 score"] >=0.5, "Pathogenic", "Benign")
        df = df.dropna()
        print(df["MutPred2 score"])
        df.to_csv("./docs/Mutpred2_Results.tsv", sep=",", index=False)

        return df
 
    
    PhD_SNPg_df = Phd_SNPg_dataframe(Phd)
    
    
    Mutpred2_df = Mutpred2_dataframe(mutpred2) 
#Beeswarm plot
    plt.figure(figsize=(8,6))
    beeswarm_p = sns.swarmplot(data=PhD_SNPg_df, x="SCORE", y="PREDICTION") 
    plt.savefig("./Figures/PhD_SNPg_Score_Beeswarm", dpi=300)   
    plt.figure(figsize=(8,6))
    beeswarm_m = sns.swarmplot(data=Mutpred2_df, x= "MutPred2 score", y="Pathogenicity")      
    plt.savefig("./Figures/MutPred2_Score_Beeswarm", dpi=300)
    
main()

