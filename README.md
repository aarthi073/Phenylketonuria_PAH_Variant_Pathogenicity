**PAH Variant Pathogenicity Analysis**

*Overview*

This project develops a bioinformatics pipeline to analyze missense variants in the PAH gene associated with phenylketonuria (PKU). The workflow retrieves clinical variant data from ClinVar, preprocesses pathogenic missense substitutions, predicts variant pathogenicity, and visualizes mutation characteristics.


*Repository Structure*
```
.
|-- Figures
|   |-- MutPred2_Score_Beeswarm.png
|   |-- MutPred2_Score_Stripplot.png
|   |-- Phd_SNPg_Alternate_Nucleotide_Stripplot.png
|   |-- Phd_SNPg_PhyloP100_Beeswarm.png
|   |-- Phd_SNPg_Score_Beeswarm.png
|   `-- SVM_Confusion_Matrix.png
|-- Peptides
|   `-- wildtype.fasta
|-- congif.yaml
|-- docs
|   |-- Mutpred2_Results.tsv
|   |-- Phd-Snp_Results.tsv
|   |-- Phd_SNPg_Results.txt
|   |-- cleaned.json
|   |-- clinvar_result.txt
|   |-- esearch.json
|   |-- esummary.json
|   |-- features.txt
|   |-- label.txt
|   |-- mutpred2_results.txt
|   |-- substitutions.txt
|   `-- test.py
|-- scripts
|   |-- SVM.py
|   |-- accuracy.py
|   |-- beeswarm.py
|   |-- boxplot.py
|   |-- clean.py
|   |-- esummary.py
|   |-- mutpred.sh
|   `-- variants_esearch.py
`-- README.md

```
*Objectives*
Retrieve PAH variants from ClinVar
Filter for missense mutations
Predict pathogenicity using MutPred2 (run using Linux skills) and PhD-SNPg (run on web server)
Visualize mutation distributions
Build machine learning models to classify pathogenic variants

*Workflow*
1. Retrieve ClinVar records through the NCBI E-utilities API
2. Parse XML records using Biopython
3. Filter missense variants
4. Extract HGVS protein substitutions
5. Run pathogenicity prediction tools
6. Generate visualizations
7. Train machine learning classifiers (in progress)

*Technologies*
Python
Biopython
pandas
requests
Conda
Linux
ClinVar
MutPred2
PhD-SNPg

*Current Status*
Work in Progress

*Upcoming additions*

Random Forest classifier
Support Vector Machine classifier
Feature engineering using amino acid physicochemical properties like isoelectric point (pI) and hydrophobicity.
Pipeline automation with Snakemake
