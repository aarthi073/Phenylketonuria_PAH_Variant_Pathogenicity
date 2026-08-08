import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def main():
    try:
        mutpred2 = sys.argv[1]
    except IndexError:
            print("no pathogenicity prediction table passed")



    df = pd.read_csv(mutpred2)      
       
    box = sns.boxplot(data=df, x="MutPred2 score", y="Pathogenicity")
    plt.show()   
main()
