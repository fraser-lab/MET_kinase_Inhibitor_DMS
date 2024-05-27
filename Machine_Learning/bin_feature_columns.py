import pandas as pd
import numpy as np
import re


inhibitors = ["Crizo", "Tepo", "NVP", "Mere", "Savo", "Cabo", "Camp", "Gle", "Glu"] # "Tiv","A458","DMSO"

def bin_dataframe_columns(original_df):
    original_df_columns = list(original_df.columns.values)
    dddG_bins = [-4,-1,1,4]
    vol_bins = [-170,-100,100,170]
    dist_bins = [5,10,30,40]
    mol_bins = [300,400,550,650]
    inhib_dist_bins = [0,5,25,35]
    rmsf_bins = [0,3,6,9]
    rfscore_bins=[5.9,6.1,6.3,6.6]
    pocket_vol_bins = [300,800,1300,1800,2600]
    hydro_bins = [-13,10,30,70]
    polar_bins = [3,11,16,24]
    rmsd_bins = [0,2,10,30]
    ligrmsd_bins = [0,6,10,20,35]
    dddG_bin_labels = [f'Bin {i+1}' for i in range(len(dddG_bins)-1)]
    vol_bin_labels = [f'Bin {i+1}' for i in range(len(vol_bins)-1)]
    dist_bin_labels = [f'Bin {i+1}' for i in range(len(dist_bins)-1)]
    inhib_dist_bin_labels = [f'Bin {i+1}' for i in range(len(inhib_dist_bins)-1)]
    mol_bin_labels = [f'Bin {i+1}' for i in range(len(mol_bins)-1)]
    rmsf_bin_labels = [f'Bin {i+1}' for i in range(len(rmsf_bins)-1)]
    rfscore_bin_labels = [f'Bin {i+1}' for i in range(len(rfscore_bins)-1)]
    pocket_vol_bin_labels = [f'Bin {i+1}' for i in range(len(pocket_vol_bins)-1)]
    hydro_bin_labels = [f'Bin {i+1}' for i in range(len(hydro_bins)-1)]
    polar_bin_labels = [f'Bin {i+1}' for i in range(len(polar_bins)-1)]
    rmsd_bin_labels = [f'Bin {i+1}' for i in range(len(rmsd_bins)-1)]
    ligrmsd_bin_labels = [f'Bin {i+1}' for i in range(len(ligrmsd_bins)-1)]
    p2rank_bin_labels = [f'Bin {i+1}' for i in range(len(p2rank_bins)-1)]
    original_df['dddG_Bin'] = pd.cut(original_df['dddG'], bins=dddG_bins, labels=dddG_bin_labels, include_lowest=True, right=False)
    original_df['dddG'] = original_df.groupby('dddG_Bin')['dddG'].transform('median')
    original_df['vol_Bin'] = pd.cut(original_df['volume_difference'], bins=vol_bins, labels=vol_bin_labels, include_lowest=True, right=False)
    original_df['volume_difference'] = original_df.groupby('vol_Bin')['volume_difference'].transform('median')
    original_df['dist_Bin'] = pd.cut(original_df['distance'], bins=dist_bins, labels=dist_bin_labels, include_lowest=True, right=False)
    original_df['distance'] = original_df.groupby('dist_Bin')['distance'].transform('median')
    original_df['inhib_dist_Bin'] = pd.cut(original_df['inhib_distance'], bins=inhib_dist_bins, labels=inhib_dist_bin_labels, include_lowest=True, right=False)
    original_df['inhib_distance'] = original_df.groupby('inhib_dist_Bin')['inhib_distance'].transform('median')
    original_df['mol_Bin'] = pd.cut(original_df['molecular_weight'], bins=mol_bins, labels=mol_bin_labels, include_lowest=True, right=False)
    original_df['molecular_weight'] = original_df.groupby('mol_Bin')['molecular_weight'].transform('median')
    original_df['rmsf_Bin'] = pd.cut(original_df['rmsf'], bins=rmsf_bins, labels=rmsf_bin_labels, include_lowest=True, right=False)
    original_df['rmsf'] = original_df.groupby('rmsf_Bin')['rmsf'].transform('median')
    original_df['rfscore_Bin'] = pd.cut(original_df['rfscore'], bins=rfscore_bins, labels=rfscore_bin_labels, include_lowest=True, right=False)
    original_df['rfscore'] = original_df.groupby('rfscore_Bin')['rfscore'].transform('median')
    original_df['pocket_Bin'] = pd.cut(original_df['pocket_volume'], bins=pocket_vol_bins, labels=pocket_vol_bin_labels, include_lowest=True, right=False)
    original_df['pocket_volume'] = original_df.groupby('pocket_Bin')['pocket_volume'].transform('median')
    original_df['hydro_Bin'] = pd.cut(original_df['hydrophobicity_score'], bins=hydro_bins, labels=hydro_bin_labels, include_lowest=True, right=False)
    original_df['hydrophobicity_score'] = original_df.groupby('hydro_Bin')['hydrophobicity_score'].transform('median')
    original_df['polar_Bin'] = pd.cut(original_df['polarity_score'], bins=polar_bins, labels=polar_bin_labels, include_lowest=True, right=False)
    original_df['polarity_score'] = original_df.groupby('polar_Bin')['polarity_score'].transform('median')
    original_df['rmsd_Bin'] = pd.cut(original_df['rmsd'], bins=rmsd_bins, labels=rmsd_bin_labels, include_lowest=True, right=False)
    original_df['rmsd'] = original_df.groupby('rmsd_Bin')['rmsd'].transform('median')
    original_df['ligrmsd_Bin'] = pd.cut(original_df['ligrmsd'], bins=ligrmsd_bins, labels=ligrmsd_bin_labels, include_lowest=True, right=False)
    original_df['ligrmsd'] = original_df.groupby('ligrmsd_Bin')['ligrmsd'].transform('median')
    original_df['p2rank_Bin'] = pd.cut(original_df['p2rank_score'], bins=p2rank_bins, labels=p2rank_bin_labels, include_lowest=True, right=False)
    original_df['p2rank_score'] = original_df.groupby('p2rank_Bin')['p2rank_score'].transform('median')

    original_df = original_df.drop(['dddG_Bin','vol_Bin','dist_Bin','mol_Bin','inhib_dist_Bin','rmsf_Bin','rfscore_Bin','pocket_Bin','hydro_Bin','polar_Bin','rmsd_Bin','ligrmsd_Bin','p2rank_Bin'],axis=1)
    original_df = original_df[original_df_columns]
    return original_df

data_df = pd.read_csv('data/all_features_all_data_pre_binning.csv',header=0)
print(data_df.columns)
data_df = bin_dataframe_columns(data_df)
data_df.to_csv("data/all_features_all_data.csv",index = False)