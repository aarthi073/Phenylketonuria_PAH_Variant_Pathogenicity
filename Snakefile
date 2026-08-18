rule all:
        input:
                "Figures/SVM_PCA.png",
                "Figures/Pathogenicity_Confusion_Matrix_SVM.png",
                "Figures/Important_Features_SVM.png",
                "Figures/RF_PCA.png",
                "Figures/Pathogenicity_Confusion_Matrix_Random_Forest.png",
                "Figures/Important_Features_Random_Forest.png"

rule clinVar_results:
        output:
                json_file = "docs/esearch.json"
        script:
                "scripts/variants_esearch.py"

rule result_summary:
        output:
                json_file = "docs/esummary.json"
        script:
                "scripts/esummary.py"
rule plots:
        input:
                phd_snp = "docs/Phd_SNPg_Results.txt",
                mutpred2 = "docs/mutpred2_results.txt"
        output:
                "Figures/MutPred2_Score_Beeswarm.png",
                "Figures/PhD_SNPg_Score_Beeswarm.png"
        script:
                "scripts/beeswarm.py"

rule SVM:
        output: 
               "Figures/SVM_PCA.png",
               "Figures/Pathogenicity_Confusion_Matrix_SVM.png",
               "Figures/Important_Features_SVM.png"
        script:
               "scripts/SVM.py"
          
rule RF:
        output:
               "Figures/RF_PCA.png",
               "Figures/Pathogenicity_Confusion_Matrix_Random_Forest.png",
               "Figures/Important_Features_Random_Forest.png"
        script:
               "scripts/RF.py"

