import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from tqdm import tqdm

import clustercv


# Read data file
data = pd.read_csv("cmpd_activity_tcrd_5.4.1.csv", low_memory=False)
# Filter data
targets = clustercv.filter_data(data, act_type='IC50', min_samples=100)
# Generate features
scores = []
for i in tqdm(range(len(targets))):
    n = targets[i]
    target_data = clustercv.get_target_data(data, n)
    features = clustercv.batchECFP(target_data.smiles)
    # Cluster features and generate folds
    n_folds = 10
    folds, max_score = clustercv.ClusterCV(features, n_folds=10)
    scores.append(max_score)
    
plt.hist(scores)
plt.xlabel("Silhouette score")
plt.ylabel("Frequency")
plt.savefig("SilhouetteScores.png", format='pdf')