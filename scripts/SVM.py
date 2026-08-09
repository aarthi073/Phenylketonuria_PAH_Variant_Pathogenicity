import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from sklearn.svm import LinearSVC,SVC
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
                else:
                
                        orig_aa.append(match.group(1))
                        pos.append(match.group(2))
                        new_aa.append(match.group(3))
        return pd.DataFrame({
        "orig_aa": orig_aa,       
        "pos":pos,
        "new_aa":new_aa 
        })          

df_new = subs(df)


counts = df["pathogenicity"]["pathogenicity"].value_counts()
unusable_classes = counts[counts < 2].index
if len(unusable_classes) > 0:
        for c in unusable_classes: 
                df["pathogenicity"]["pathogenicity"] = df["pathogenicity"]["pathogenicity"].replace(f"{c}", "Benign/Likely benign")


#define features (x) and response variable (y)
X = df_new[["orig_aa", "pos", "new_aa"]]
y = df["pathogenicity"]

#Label Encoding: transforms text categories of the response variable into numbers
le = LabelEncoder()
#1D array expected
y_encoded = le.fit_transform(np.ravel(y))

#print(unusable_classes)
#print(le.classes_)
#print(pd.Series(y_encoded).value_counts())	

#One-hot encoding and standardizing features
#multi-columns of 0/1 for each amino acid 

preprocessor = ColumnTransformer(
        transformers = [
            ("cat", OneHotEncoder(handle_unknown="ignore"), ["orig_aa", "new_aa","pos"])
            #("num", StandardScaler(), ["pos"])
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
#plt.show()



#train-test split
X_temp, X_test, y_temp, y_test = train_test_split(X_processed, y_encoded, test_size=0.15, stratify=y_encoded,random_state=42)
X_train, X_val, y_train, y_val = train_test_split(X_temp, y_temp, test_size=0.15, stratify=y_temp, random_state=42)
#Train a SVM classifier
# Standard tuning grid for an RBF SVM
# RBF-SVM: Accuracy Evaluation Metric
C_values = [0.1, 1, 10, 100, 1000]
gamma_values = [0.01, 0.1, 0.5, 1, 5]

best_score = 0
best_train_score = 0
best_params = None
best_train_params = None

# HPO (grid search method)
for C in C_values:
    for gamma in gamma_values:
        model = SVC(kernel="rbf", C=C, gamma=gamma)
        model.fit(X_train, y_train)
        
        y_pred_val = model.predict(X_val)
        val_acc = accuracy_score(y_val, y_pred_val)
        print(f"C={C}, gamma={gamma}, val_acc={val_acc:.3f}")

         # training accuracy
        y_pred_train = model.predict(X_train)
        train_acc = accuracy_score(y_train, y_pred_train)
        print(f"C={C}, gamma={gamma}",
              f"train_acc={train_acc:.3f}")
        
        if val_acc > best_score:
            # update the score (best_score)
            best_score = val_acc
            # update the hyperparameters (best_params)
            best_params = (C, gamma)
            #update the training accuracy associated with the best validation score
            best_val_train_acc = train_acc
        if train_acc > best_train_score:
            # update the score (best_train_score)
            best_train_score = train_acc
            # update the hyperparameters (best_train_params)
            best_train_params = (C, gamma)

# --- 4. Display the best parameters ---
print("\n Best parameters:")
print(f"C = {best_params[0]}, gamma = {best_params[1]}, validation accuracy = {best_score:.3f}, training accuracy = {best_val_train_acc:.3f}")
svm_acc = best_score

#print(f"C = {best_train_params[0]}, gamma = {best_train_params[1]}, training accuracy = {best_train_score:.3f}")
svm_acc = best_score


#svm_poly = SVC(kernel="poly", class_weight="balanced",degree=3, coef0=1, C=1)
#svm_poly.fit(X_train, y_train)
#y_pred = svm_poly.predict(X_test)


#Accuracy
#acc = accuracy_score(y_test, y_pred, normalize=True, sample_weight=None)
#print(f"accuracy: {acc:.3f}")


#Final Model Based on Best Accuracy
final_model = SVC(kernel="rbf", C=best_params[0], gamma = best_params[1]).fit(X_train, y_train)
y_pred_test = final_model.predict(X_test)
test_acc = accuracy_score(y_test, y_pred_test)
print(f"Final Test Accuracy (svm_rbf): {test_acc:.3f}")
#visualizations: confusion matrix heatmaps, amino acid substitution frequency heatmaps.
cm = confusion_matrix(y_test, y_pred_test)
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Pathogenicity Confusion Matrix: SVM")
plt.show()


#Amino Acid Substitution Heatmap
#aa_df = pd.DataFrame(df, 
#plt.figure(figsize=(10,8))


