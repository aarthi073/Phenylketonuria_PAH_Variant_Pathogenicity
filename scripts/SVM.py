import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.svm import LinearSVC,SVC
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix
from sklearn.compose import ColumnTransformer
import pandas as pd
import os 

docs = "../docs"

files = { 
         "substitutions": f"{docs}/substitutions.txt",
         "pathogenicity": "./label.txt"
}

df = {}
for key, path_value in files.items():

  with open(path_value, "r") as f:
                lines = [line.strip() for line in f]
                df1 = pd.DataFrame(lines, columns=[key])
  df[key] = df1


def subs(X):
        orig_aa = []
        pos = []
        new_aa = []
        for i in X["substitutions"]["substitutions"]:
                match = re.match(r"([A-Z])(\d+)([A-Z])",i)
                if not match:
                        print(f"Failed to parse: {i}")
                orig_aa.append(match.group(1))
                pos.append(match.group(2))
                new_aa.append(match.group(3))
        return pd.DataFrame({
        "orig_aa": orig_aa,       
        "pos":pos,
        "new_aa":new_aa 
        })          

df_new = subs(df)

#define features (x) and response variable (y)
X = df_new[["orig_aa", "pos", "new_aa"]]
y = df["pathogenicity"]

#Label Encoding: transforms text categories of the response variable into numbers 
le = LabelEncoder()
#1D array expected
y_encoded = le.fit_transform(np.ravel(y))

print(le.classes_)
#One-hot encoding and standardizing features
#multi-columns of 0/1 for each amino acid 

preprocessor = ColumnTransformer(
        transformers = [
            ("cat", OneHotEncoder(handle_unknown="ignore"), ["orig_aa", "new_aa"]),
            ("num", StandardScaler(), ["pos"])
        ]
)

X_processed = preprocessor.fit_transform(X)

#train-test split
X_train, X_test, y_train, y_test = train_test_split(X_processed, y_encoded, test_size=0.2, random_state=42)

#Train a SVM classifier
svm_poly = SVC(kernel="poly", class_weight="balanced",degree=3, coef0=1, C=1)
svm_poly.fit(X_train, y_train)
y_pred = svm_poly.predict(X_test)

#visualizations: confusion matrix heatmaps, amino acid substitution frequency heatmaps.
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix")
plt.show()
