import numpy as np
import pandas as pd
import sys, os, pickle
import argparse, textwrap
from preprocess import data_preparation
from common import  build_in_vivo_model, transfer_test_on_in_vitro

def main():
    
    parser = argparse.ArgumentParser(description = 'Pipeline for buidling Malaria Artemisinin resistance in vivo prediction models and transfer learning on in vitro datasets.',
            usage = 'use "python %(prog)s --help" for more information',
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--train_path', type = str, 
            help = "path to your training data, in .csv format", 
            default = "../../rawdata/SubCh2_TrainingData.csv")
    parser.add_argument('--valid_path', type = str, 
            help = "path to your transfer validation data, in .csv format", 
            default = "../../rawdata/SubCh1_TrainingData.csv")
    parser.add_argument('-m','--model_type', type = str, 
            help = '''machine learning models to use:
            lgb: LightGBM; 
            xgb: XGBoost;
            rf: Random Forest; 
            gpr: Gaussian Process Regression;
            lr: Linear Regression;
            default: lgb''', 
            default = 'lgb')
    parser.add_argument('--no_quantile',
            action = "store_true",
            help = "if specified, do not use quantile normalization.")
    parser.add_argument('-n','--use_top_features', type = int, 
            help = '''if specified, used top features based on shap analysis.
            for example: 3, 5, 10, 20 or 30; 
            Otherwise, use the whole genome''',
            default = None)
    parser.add_argument('--feature_selection_mode', 
            action = "store_true",
            help = "if used, begins leave-one-out feature selection strategy from specified gene set")

    args = parser.parse_args()
    
    if args.feature_selection_mode:
        print("Begin leave-one-out feature selection mode ...")
        assert args.use_top_features != None, "Must specify -n for leave-one-out feature selection mode!"
    if args.no_quantile:
        print("Do not use normalization ...")
    opts = vars(parser.parse_args())
    
    run(**opts)

def leave_one_out_selection(n):
    """ Generate features to use when leave-one-out feature selection strategy is on
    
    Params
    -------
    n: int
        number of top genes used in feature selection

    Yields
    -------
    common_genes: a list of str
        str must be present in training and testing feature sets

    """
    path = '../downstream_analysis/invivo_top_'+str(n)+'_genes.pkl'
    genes = pickle.load(open(path, 'rb'))
    for g in genes:
        remain = genes.copy()
        remain.remove(g)
        print("Leave-one-out:", g, "Remaining:", remain)
        yield(g, remain)


def run(train_path, valid_path, model_type, no_quantile, use_top_features, feature_selection_mode):
    
    # load training dataset
    df_train = pd.read_csv(train_path) # in vivo dataset for building the model
    # load transfering test dataset
    df_val = pd.read_csv(valid_path) # in vitro dataset for transfer testing
    
    # preprocess both in vivo and in vitro dataset
    if feature_selection_mode:
        for g, genes in leave_one_out_selection(use_top_features):
            df_invivo, df_invitro = data_preparation(df_train, df_val, genes, feature_selection_mode)
            generate_results(df_invivo, df_invitro, model_type, g+'_')
    else:
        if use_top_features != None:
            path = '../downstream_analysis/invivo_top_'+str(use_top_features)+'_genes.pkl'
            genes = pickle.load(open(path, 'rb'))
        else:
            genes = None
        df_invivo, df_invitro = data_preparation(df_train, df_val, no_quantile, genes)
        generate_results(df_invivo, df_invitro, model_type, '')

def generate_results(df_invivo, df_invitro, model_type, path_affix):
    """
    Params
    ------
    df_invivo: a Pandas DataFrame
    df_invitro: a Pandas DataFrame
    model_type: str
    path_affix: str
    """
    #Part 1: five fold cross validation on in vivo data
    eva_cv, eva_cv_conf = build_in_vivo_model(df_invivo, model_type)
    # Part 2: transfer testing on in vitro data
    eva_tv, eva_tv_conf = transfer_test_on_in_vitro(df_invitro)

    # save evaluation results
    path = './'+path_affix+'performance'
    os.makedirs(path, exist_ok = True)
    eva_cv.to_csv(path+'/in_vivo_cv_results.csv', index = False)
    eva_cv_conf.to_csv(path+'/in_vivo_cv_confidence.csv', index = False)
    eva_tv.to_csv(path+'/in_vitro_tv_results.csv', index = False)
    eva_tv_conf.to_csv(path+'/in_vitro_tv_confidence.csv', index = False)

if __name__ == "__main__":
    main()
