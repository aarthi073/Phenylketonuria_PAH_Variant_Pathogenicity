import matplotlib as plt
def main():
    try:
        mutpred2 = sys.argv[1]
    except IndexError:
            print("no pathogenicity prediction table passed")


    def Phd_SNPg_dataframe(Phd_SNPg):
        headers = ["CHROM", "REF", "ALT", "MUT", "CODING", "PREDICTION", "SCORE", "PhyloP100"]
        cols = ["CHROM", "POS", "ID", "REF", "ALT", "MUT", "CODING", "PREDICTION", "SCORE", "FDR", "Ph>
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
        df = df.drop(columns=["POS", "ID", "FDR", "AvgPhyloP100",  "transcript", "gene", "strand", "co>
        df.to_csv("../docs/Phd-Snp_Results.tsv", sep="\t", index=False)
        return df

