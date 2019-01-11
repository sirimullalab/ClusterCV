from random import shuffle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from sklearn.neighbors import DistanceMetric
from sklearn.metrics import silhouette_score
from sklearn.cluster import AgglomerativeClustering


def filter_data(data, act_type='IC50', min_samples=200):
    """Returns a list of target_ids that meet the given criteria."""
    # Remove items with missing smiles
    data = data[[not f for f in data.smiles.isna()]]
    # Filter by activity type
    if act_type is not None:
        data = data[data.act_type == act_type]
    # Filter by number of samples per target
    targets = np.unique(data.target_id)
    target_ids = []
    for t in targets:
        if (data[data.target_id == t].shape[0] >= min_samples):
            target_ids.append(t)
    return target_ids


def get_target_data(data, target_id, act_type='IC50'):
    """Returns a data frame of all the ligands for a given target
       Also makes sure that all the smiles are valid, and
       filters by weight.""" 
    if act_type is not None:
        data = data[data.act_type == act_type]
    target_data = data[data.target_id == target_id]
    # Filter by molecules that can be converted by rdkit
    n_ligs =  target_data.shape[0]
    mols = np.zeros(n_ligs, dtype=object)
    for i in range(n_ligs):
        try:
            mols[i] = MolFromSmiles(target_data.smiles.iloc[i])
        except:
            mols[i] = None
    mols = pd.Series(mols)
    target_data = target_data[[not m for m in mols.isna()]]
    # Filter by weight
    weights = target_data.smiles.apply(lambda x: ExactMolWt(MolFromSmiles(x)))
    target_data = target_data[(weights >= 100) & (weights <= 600)]
    return target_data


def batchECFP(smiles, radius=4, nBits=1024):
    """Calculate ECFP fingerprints for a list of smiles.
    Parameters:
    smiles - An iterable with SMILES strings.
    radius - the fingerprint radius
    nBits - the number of bits to use to encode
    Returns: A feature matrix.
    """
    smiles = np.array(smiles)
    n = len(smiles)
    fingerprints = np.zeros((n, nBits), dtype=int)
    for i in range(n):
        mol = MolFromSmiles(smiles[i])
        fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        fingerprints[i] = np.array(list(fp.ToBitString()))
    return(fingerprints)


def ClusterCV(features, linkage='single', n_folds=10):
    """ Cluster-cross-validation.
    linkage - 'single', 'complete', or 'average'
    """
    # Get distance matrix
    dist = DistanceMetric.get_metric('jaccard')
    distmatrix = dist.pairwise(features)
    # Find best number of clusters
    num_samples = features.shape[0]
    scores = np.zeros(num_samples)
    for i in range(n_folds * 2, num_samples):
        clustering = AgglomerativeClustering(n_clusters=i,
                                             affinity='precomputed',
                                             linkage=linkage)
        clustering = clustering.fit(distmatrix)
        scores[i] = silhouette_score(distmatrix, labels=clustering.labels_,
                                     metric='precomputed')
    max_score = max(scores)
    n_clusters = np.where(scores == max(scores))[0][0]
    #print(scores)
    print("Number of clusters:", n_clusters)
    # Cluster
    clustering = AgglomerativeClustering(n_clusters=n_clusters,
                                         affinity='precomputed',
                                         linkage=linkage)
    clustering = clustering.fit(distmatrix)
    # Randomly assign each cluster to one fold
    cluster_nums = list(range(n_clusters))
    shuffle(cluster_nums)
    fold_assignments = {}
    for n in range(n_clusters):
        fold_assignments[cluster_nums[n]] = n % n_folds
    # Assign samples to folds
    folds = np.zeros(num_samples)
    for j in range(num_samples):
        folds[j] = fold_assignments[clustering.labels_[j]]
    return folds, max_score


if __name__ == "__main__":
    import sys

    data = pd.read_csv("cmpd_activity_tcrd_5.4.1.csv", low_memory=False)
    targets = filter_data(data)
    n = targets[int(sys.argv[1])]
    target_data = get_target_data(data, n)
    print(target_data.smiles)
    features = batchECFP(target_data.smiles)
    plt.matshow(features)
    plt.show()
    folds = ClusterCV(features)
    print(folds)
