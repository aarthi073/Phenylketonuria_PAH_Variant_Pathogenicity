import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.compose import ColumnTransformer
import pandas as pd
import os



docs = "../docs"

files = {
         "substitutions": f"{docs}/substitutions.txt",
         "pathogenicity": f"{docs}/label.txt"
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
            ("cat", OneHotEncoder(handle_unknown="ignore"), ["orig_aa", "new_aa","pos"])
            #("num", StandardScaler(), ["pos"])
        ]
)

X_processed = preprocessor.fit_transform(X)

#train-test, validation splits
X_temp, X_test, y_temp, y_test = train_test_split(X_processed, y_encoded, test_size=0.15, stratify=y_encoded, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.15, stratify=y_temp, random_state=42)
#Random Forest model training: n_estimators: 5, 10, 50, 100, 500
candidates = [5, 10, 50, 100, 500]
results = []

for n in candidates:
        RF = RandomForestClassifier(n_estimators=n, max_depth=None, random_state=42)
        RF.fit(X_train, y_train)
        y_pred_val = RF.predict(X_val)
        acc= accuracy_score(y_pred_val, y_val)
        results.append(acc)
        print(f"estimators={n} accuracy: {acc:.3f}")
best_n = candidates[np.argmax(results)] 
print(f"Best number of estimators by validation: {best_n}")


#Visualizations
#Confusion Matrix
cm = confusion_matrix(y_val, y_pred_val)
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Pathogenicity Confusion Matrix: Random Forest")
plt.show()
