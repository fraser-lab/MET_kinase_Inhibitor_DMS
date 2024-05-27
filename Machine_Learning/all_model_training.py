#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:41:28 2024

@author: ashraya
"""

from itertools import combinations
import xgboost as xgb
import re
import sys
import numpy as np
import pandas as pd
import pickle as pkl
from scipy import stats
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from warnings import filterwarnings
filterwarnings('ignore')
base_dir = "/wynton/group/fraser/aravikumar/MET_kinase_Inhibitor_DMS/Machine_Learning/"

feature_list = {"esm":['score','1','ESM LLR'],
                "dddg":['dddG','0','∆∆∆G'],
                "ddg_all":['ddG_all','-1','∆∆G'],
                "atp_distance":['distance','0','Residue->ATP Distance'],
                "volume_difference":['volume_difference','0','∆Volume'],
                "inhibitor_weight":['molecular_weight','0','MWt'],
                "inhibitor_distance":['inhib_distance','0','Inhibitor distance'],
                "rf_score":['rfscore','0','RF Score'],
                "crystal_rmsf":['rmsf','0','MET Crystal RMSF'],
                "pocket_volume":['pocket_volume','0','Pocket Volume'],
                "hydrophobicity":['hydrophobicity_score','0','Hydrophobicity Score'],
                "polarity":['polarity_score','0','Polarity Score'],
                "residue_rmsd":['rmsd','0','Residue RMSD'],
                "ligand_rmsd":['ligrmsd','0','Ligand RMSD'],    
    }
hyperparam_dict=dict()
list_of_models =[]
infile=open(base_dir+"data/best_hyperparameters_final.csv")
for line in infile:
    lineparts=line.split("|")
    hyp_dict = eval(lineparts[1])
    list_of_models.append(lineparts[0])
    hyperparam_dict[lineparts[0]]=hyp_dict
infile.close()

def generate_combinations(arr):
    all_combinations = []
    n = len(arr)
    for r in range(1, n + 1):
        combinations_of_length_r = combinations(arr, r)
        all_combinations.extend(combinations_of_length_r)
    return all_combinations

def prepare_model_dict(features):
    model_dict=dict()
    monotone_array=[]
    feat_array=[]
    label_array=[]
    for feature in features:
        column_name,monotone,label = feature_list[feature]
        label_array.append(feature)
        monotone_array.append(monotone)
        feat_array.append(column_name)
    monotone_string="("
    monotone_string+=",".join(x for x in monotone_array)
    monotone_string+=")"
    dict_key="_".join(x for x in label_array)
    model_dict[dict_key]=[feat_array,monotone_string]
    return dict_key,model_dict
  

def prepare_feature_string(features):
    model_label=[]
    for feat in features:
        model_label.append(feature_list[feat][2])
    return "+".join(x for x in model_label)


def get_holdout_df(combined_df):
    keys = ["Crizo", "Gle"]
    holdout_df = combined_df[combined_df['key'].isin(keys)].copy()
    train_df = combined_df[~combined_df['key'].isin(keys)].copy()

    chooser = np.random.default_rng(7)
    holdout_positions = chooser.choice(train_df['pos'].unique(), size=int(len(train_df['pos'].unique()) * 0.2), replace=False)
    holdout_positions_df = train_df[train_df['pos'].isin(holdout_positions)]
    train_df = train_df[~train_df['pos'].isin(holdout_positions)]

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    holdout_mutations = chooser.choice(amino_acids, size=int(len(amino_acids) * 0.1), replace=False)
    train_df['mut'] = train_df['pos_mut'].apply(lambda x: re.sub(r'\d+(\w)', r'\1', x))
    holdout_mutations_df = train_df[train_df['mut'].isin(holdout_mutations)].drop(['mut'], axis=1)
    train_df = train_df[~train_df['mut'].isin(holdout_mutations)].drop(['mut'], axis=1)

    holdout_df = pd.concat([holdout_positions_df, holdout_mutations_df, holdout_df])
    return train_df,holdout_df

def cross_validate_model(model_name,model_dict,full_train_df,hyperparameters):
    train_on = model_dict[model_name][0]
    hyperparameters['monotone_constraints'] = model_dict[model_name][1]
        
    to_drop = [column for column in full_train_df.columns if column not in train_on]
    model = xgb.XGBRegressor(**hyperparameters)
    
    mse_train_scores = []
    pearson_train_scores = []
    mse_validate_scores = []
    pearson_validate_scores = []
    
    all_positions = full_train_df['pos'].unique()
    full_train_df['mut'] = full_train_df['pos_mut'].apply(lambda x: re.sub(r'\d+(\w)', r'\1', x))
    all_aa = full_train_df['mut'].unique()
    full_train_df = full_train_df.drop(['mut'], axis = 1)
    rng = np.random.default_rng(120)
    rng.shuffle(all_positions)
    rng.shuffle(all_aa)
    split_positions = np.split(all_positions, 10)
    split_aa = np.split(all_aa,9)
    
    #Cross validation loop
    for i in range(10):
        temp_train_df = full_train_df.copy()
        holdout_positions = split_positions[i]
        holdout_positions_df = temp_train_df[temp_train_df['pos'].isin(holdout_positions)]
        temp_train_df = temp_train_df[~temp_train_df['pos'].isin(holdout_positions)]
        if i > 8:
            aa_index = rng.choice(range(8),size =1)[0]
        else:
            aa_index = i
        holdout_mutations = split_aa[aa_index]
        temp_train_df['mut'] = temp_train_df['pos_mut'].apply(lambda x: re.sub(r'\d+(\w)', r'\1', x))
        holdout_mutations_df = temp_train_df[temp_train_df['mut'].isin(holdout_mutations)].drop(['mut'], axis=1)
        temp_train_df = temp_train_df[~temp_train_df['mut'].isin(holdout_mutations)].drop(['mut'], axis=1)
        validate_df = pd.concat([holdout_positions_df, holdout_mutations_df])
        

        X_train = temp_train_df.drop(to_drop, axis=1)
        X_train = X_train[train_on]
        y_train = temp_train_df['mean']
        X_test = validate_df.drop(to_drop, axis=1)
        X_test = X_test[train_on]
        y_test = validate_df['mean']

        # Train
        model.fit(X_train, y_train)
        
        y_pred = model.predict(X_test)
        y_train_pred = model.predict(X_train)
        mse_train_scores.append(mean_squared_error(y_train, y_train_pred))
        pearson_train_scores.append(stats.pearsonr(y_train, y_train_pred)[0])
        mse_validate_scores.append(mean_squared_error(y_test, y_pred))
        pearson_validate_scores.append(stats.pearsonr(y_test, y_pred)[0])
    return mse_train_scores, pearson_train_scores, mse_validate_scores, pearson_validate_scores

def train_and_test_model(model_name,model_dict,train_df,test_df,hyperparameters):
    
    train_on = model_dict[model_name][0]
    hyperparameters['monotone_constraints'] = model_dict[model_name][1]
    to_drop = [column for column in train_df.columns if column not in train_on]    
    model = xgb.XGBRegressor(**hyperparameters)
    X_train = train_df.drop(to_drop, axis=1)
    X_train = X_train[train_on]
    y_train = train_df['mean']
    X_test = test_df.drop(to_drop, axis=1)
    X_test = X_test[train_on]
    y_test = test_df['mean']
    print(X_train.columns)
    print(hyperparameters['monotone_constraints'])
    
    # Train
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    y_train_pred = model.predict(X_train)
    mse_train_score = mean_squared_error(y_train, y_train_pred)
    pearson_train_score = stats.pearsonr(y_train, y_train_pred)[0]
    mse_test_score = mean_squared_error(y_test, y_pred)
    pearson_test_score = stats.pearsonr(y_test, y_pred)[0]
    return model,mse_train_score, pearson_train_score, mse_test_score, pearson_test_score

data_df = pd.read_csv(base_dir+"data/all_features_all_data.csv",header = 0)
train_df,test_df = get_holdout_df(data_df)
feature_set=sys.argv[1].split(",")
print(feature_set)
model_name,model_dict=prepare_model_dict(feature_set)
best_parameters = hyperparam_dict[model_name]
mse_train_scores, pearson_train_scores, mse_validate_scores, pearson_validate_scores = cross_validate_model(model_name,model_dict, train_df, best_parameters)
print("Cross validation results:"+prepare_feature_string(feature_set)+","+str(np.average(mse_train_scores))+","+str(np.average(mse_validate_scores))+","+str(np.average(pearson_train_scores))+","+str(np.average(pearson_validate_scores))+","+str(np.std(pearson_train_scores))+","+str(np.std(mse_validate_scores))+","+str(np.std(pearson_validate_scores))+","+str(np.std(mse_train_scores))+"\n")
train_df = train_df.drop(['mut'],axis=1)
trained_model, mse_train_score, pearson_train_score, mse_test_score, pearson_test_score =  train_and_test_model(model_name,model_dict,train_df,test_df,best_parameters)
print("Train and test results:"+prepare_feature_string(feature_set)+","+str(mse_train_score)+","+str(pearson_train_score)+","+str(mse_test_score)+","+str(pearson_test_score)+"\n")
pkl.dump(trained_model, open(base_dir+"xgb_models_new/"+model_name+".pkl", "wb"))
