#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:53:51 2024

@author: ashraya
"""

import pandas as pd
import xgboost as xgb
import optuna
import sys
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from warnings import filterwarnings
filterwarnings('ignore')

from itertools import combinations

optuna.logging.set_verbosity(optuna.logging.FATAL)

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

def objective(trial):
    
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_val, label=y_val)

    param = {
        "silent": 1,
        "monotone_constraints": monotone,
        "objective": "reg:squarederror",
        "eval_metric": "rmse",
        "booster": trial.suggest_categorical("booster", ["gbtree", "dart"]),
        "lambda": trial.suggest_loguniform("lambda", 1e-8, 1.0),
        "alpha": trial.suggest_loguniform("alpha", 1e-8, 1.0),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
        'n_estimators': trial.suggest_categorical('n_estimators', [50, 100, 150, 200]),
    }

    if param["booster"] == "gbtree" or param["booster"] == "dart":
        param["max_depth"] = trial.suggest_int("max_depth", 1, 50)
        param["eta"] = trial.suggest_loguniform("eta", 1e-8, 1.0)
        param["gamma"] = trial.suggest_loguniform("gamma", 1e-8, 1.0)
        param["grow_policy"] = trial.suggest_categorical("grow_policy", ["depthwise", "lossguide"])
    if param["booster"] == "dart":
        param["sample_type"] = trial.suggest_categorical("sample_type", ["uniform", "weighted"])
        param["normalize_type"] = trial.suggest_categorical("normalize_type", ["tree", "forest"])
        param["rate_drop"] = trial.suggest_loguniform("rate_drop", 1e-8, 1.0)
        param["skip_drop"] = trial.suggest_loguniform("skip_drop", 1e-8, 1.0)

    # Add a callback for pruning.
    pruning_callback = optuna.integration.XGBoostPruningCallback(trial, "validation-rmse")
    bst = xgb.train(param, dtrain, evals=[(dtest, "validation")], callbacks=[pruning_callback])
    preds = bst.predict(dtest)
    # pred_labels = np.rint(preds)
    mse = mean_squared_error(y_val, preds)
    return mse

if __name__ == "__main__":
    
    feature_set=sys.argv[1].split(",")
    print(feature_set)
    train_df = pd.read_csv(base_dir+"data/full_training_data.csv", header = 0)
    model_name,model_dict = prepare_model_dict(feature_set)
    train_on = model_dict[model_name][0]
    monotone = model_dict[model_name][1]
    
    to_drop = [column for column in train_df.columns if column not in train_on]
    X_train = train_df.drop(to_drop, axis=1)
    y_train = train_df['mean']

    # Validation set
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.2, random_state=7)



    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=500, show_progress_bar=True)
    print("my_result|"+model_name+"|"+str(study.best_trial.params))
