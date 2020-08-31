import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
import sys, os, pickle
from tqdm import tqdm
from preprocess import data_preparation
from common import  build_in_vivo_model, transfer_test_on_in_vitro
from statistics import compute_auroc, c_index, pearsonr_cor, boostrapping_confidence_interval

def main():
    model_type = sys.argv[1] #argument
    
    # load training dataset
    df_invivo = pd.read_csv("../../rawdata/SubCh2_TrainingData.csv") # in vivo dataset for building the model
    # load transfering test dataset
    df_invitro = pd.read_csv("../../rawdata/SubCh1_TrainingData.csv") # in vitro dataset for transfer testing
    
    # preprocess both in vivo and in vitro dataset
    df_invivo, df_invitro = data_preparation(df_invivo, df_invitro)
    
    #Part 1: five fold cross validation on in vivo data
    eva_cv = build_in_vivo_model(df_invivo, model_type)
    # Part 2: transfer testing on in vitro data
    eva_tv, eva_tv_conf = transfer_test_on_in_vitro(df_invitro)

    # save evaluation results
    path = './performance'
    os.makedirs(path, exist_ok = True)
    eva_cv.to_csv(path+'/in_vivo_cv_results.csv', index = False)
    eva_tv.to_csv(path+'/in_vitro_tv_results.csv', index = False)
    eva_tv_conf.to_csv(path+'/in_vitro_tv_confidence.csv', index = False)

if __name__ == "__main__":
    main()
