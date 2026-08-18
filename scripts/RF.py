import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.preprocessing import StandardScaler, LabelEncoder, OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.compose import ColumnTransformer
import pandas as pd
import os

files = { "substitutions": 'docs/substitutions.txt', "pathogenicity": 'docs/label.txt', "hydrophobicity":  'docs/hydrophobicity.txt'}

df = {}
for key, path_value in files.items():

  with open(path_value, "r") as f:
                lines = [line.strip() for line in f]
                if key == "substitutions":
                        orig_aa = []
                        pos = []
                        new_aa = []
                        for line in lines:
                                match = re.match(r"([A-Z])(\d+)([A-Z])", line)
                                if not match:
                                        print(f"Failed to parse: {line}")
                                else:
                                        orig_aa.append(match.group(1))
                                        pos.append(match.group(2))
                                        new_aa.append(match.group(3))

                df[key] = lines
                df["orig_aa"] = orig_aa
                df["pos"] = pos
                df["new_aa"] = new_aa
df = pd.DataFrame(df)

#convert "pos" for position to integer so that numeric feature engineering can be applied
df["pos"] = df["pos"].astype(int)


counts = df["pathogenicity"].value_counts()
unusable_classes = counts[counts < 2].index
if len(unusable_classes) > 0:
        for c in unusable_classes:
                df["pathogenicity"] = df["pathogenicity"].replace(f"{c}", "Benign/Likely benign")


#define features (x) and response variable (y)
X = df[["orig_aa", "pos", "new_aa", "hydrophobicity"]]
y = df["pathogenicity"]

#Label Encoding: transforms text categories of the response variable into numbers
le = LabelEncoder()
#1D array expected
y_encoded = le.fit_transform(np.ravel(y))


#One-hot encoding and standardizing features
#multi-columns of 0/1 for each amino acid

preprocessor = ColumnTransformer(
        transformers = [
            ("cat", OneHotEncoder(handle_unknown="ignore"), ["orig_aa", "new_aa", "hydrophobicity"]),
            ("num", StandardScaler(), ["pos"])
        ]
)

X_processed = preprocessor.fit_transform(X)

#PCA
from sklearn.decomposition import PCA

# 1. Initialize PCA to reduce data to 2 dimensions
# If X_processed is a sparse matrix (common with OneHotEncoder), convert to dense array
X_dense = X_processed.toarray() if hasattr(X_processed, "toarray") else X_processed
pca = PCA(n_components=2, random_state=42)

# 2. Fit and transform the processed feature matrix
X_pca = pca.fit_transform(X_dense)

# 3. Create a clean DataFrame for visualization
pca_df = pd.DataFrame(data=X_pca, columns=['PC1', 'PC2'])
pca_df['Class'] = le.inverse_transform(y_encoded)  # Decode numbers back to text labels

# 4. Calculate Explained Variance to see how much information we kept
var_pc1 = pca.explained_variance_ratio_[0] * 100
var_pc2 = pca.explained_variance_ratio_[1] * 100

# 5. Render the PCA Scatter Plot
plt.figure(figsize=(10, 7))
sns.scatterplot(
    data=pca_df,
    x='PC1',
    y='PC2',
    hue='Class',
    palette='Set1',
    alpha=0.7,
    edgecolor='k',
    s=50
)

# Design elements
plt.title('PCA: Amino Acid Feature Separation', fontsize=14, pad=15)
plt.xlabel(f'Principal Component 1 ({var_pc1:.1f}% Variance)', fontsize=11)
plt.ylabel(f'Principal Component 2 ({var_pc2:.1f}% Variance)', fontsize=11)
plt.legend(title='Target Class', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, linestyle=':', alpha=0.6)
plt.tight_layout()
plt.savefig("./Figures/RF_PCA", dpi=300)



X_processed = preprocessor.fit_transform(X)

#train-test, validation splits
X_temp, X_test, y_temp, y_test = train_test_split(X_processed, y_encoded, test_size=0.15, stratify=y_encoded, random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.15, stratify=y_temp, random_state=42)
#Random Forest model training: n_estimators: 5, 10, 50, 100, 500
candidates = [5, 10, 50, 100, 500]
results = []

for n in candidates:
        RF = RFC(n_estimators=n, max_depth=None, random_state=42)
        RF.fit(X_train, y_train)
        y_pred_val = RF.predict(X_val)
        acc = accuracy_score(y_pred_val, y_val)
        results.append(acc)
        print(f"estimators={n} accuracy: {acc:.3f}")
best_n = candidates[np.argmax(results)] 
print(f"Best number of estimators by validation: {best_n}")

#Final Model Based on Best Accuracy
final_model = RFC(n_estimators=best_n, max_depth=None, random_state=42).fit(X_train, y_train)
y_pred_test = final_model.predict(X_test)
test_acc = accuracy_score(y_test, y_pred_test)
print(f"Final Test Accuracy Random Forest: {test_acc:.3f}")

#Visualizations
#Confusion Matrix
cm = confusion_matrix(y_test, y_pred_test)
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Pathogenicity Confusion Matrix: Random Forest")
plt.savefig("Figures/Pathogenicity_Confusion_Matrix_Random_Forest", dpi=300)


#feature importance
plt.figure(figsize=(7,5))
features = preprocessor.get_feature_names_out() 
importances = final_model.feature_importances_ 
indices = np.argsort(importances)[::-1][:10] 
plt.barh(np.array(features)[indices][::-1], importances[indices][::-1]) 
plt.title("Important Features - Random Forest Classifier") 
plt.savefig("Figures/Important_Features_Random_Forest", dpi=300)
