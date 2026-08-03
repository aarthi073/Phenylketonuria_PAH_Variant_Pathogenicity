import csv
import pandas as pd
results = []
with open("mutpred2_results.txt") as f:
    reader = csv.reader(f)
    header = next(reader)  # skip header row
    for row in reader:
        if len(row) < 3:
            continue
        try:
            score = float(row[2])
        except ValueError:
            continue  # skips fragment/continuation lines
        results.append((row[1], score))

df = pd.DataFrame(results, columns=["Substitution", "MutPred2 score"])
print(df)
