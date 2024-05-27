#This folder contains the machine learning source code and figure generation scripts 

1. bin_feature_columns.py - This takes all_features_all_data_pre_binning.csv as input, bins the features and outputs all_features_all_data.csv file, which is the input data for cross-validation, training and testing of XGBoost models
2. hyperparameter_tuning.py - Given a list of features as input, this script does XGBoost hyperparameter tuning and outputs the best set of hyperparameters
3. all_model_training.py - Given a list of features as input, this script does cross-validation, training and testing on data in all_features_all_data.csv for that model. It reads the hyperparameters for the model from the file best_hyperparameters_final.csv in data folder.
4. ML_figures.R - This script runs on output data given in outputs folder to generate the figure panels in Figure 7 in the paper.

The complete list of models that were tuned and trained are given in list_of_models.csv. Each line is a comma separated list of features that are incorporated in each model.

Running the python scripts require the following python libraries:
1. xgboost
2. scipy
3. scikit-sklearn
4. optuna
