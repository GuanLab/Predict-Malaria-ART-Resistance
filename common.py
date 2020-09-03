import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
import os, pickle
from scipy.stats import pearsonr, spearmanr
from models import *
from statistics import *

def build_in_vivo_model(df_invivo, model_type):
    """ Build in vivo clearance rate prediction model and five-fold cross validation
    Model parameters for each fold are saved as fold_*_model.sav
    
    Parameters:
    -----------
    df_invivo: a Pandas dataframe
        the processed in vivo dataset
    model_type: str
        type of model to train
        'lgb': Light GBM model
        'xgb': XGboost model
        'rf': Random Forest
        'lr': Linear regression
        'gpr' gaussian process regression
    
    Yields:
    -------
    eva_df: a Pandas dataframe
        evaluations of in vivo models during five-fold cross validation
    """
    os.makedirs('./invivo/', exist_ok = True) # where to save the processed in vivo datasets during five-fold cross validation
    os.makedirs('./params/', exist_ok = True) # where to save the model params
    kf = KFold(n_splits=10, shuffle = True, random_state = 0)
    eva_df = {'fold':[], 'AUROC':[], 'AUPRC':[], 'Pearsonr':[], 'Spearmanr':[]}
    eva_conf_df = {'Pearsonr mean[95CI]':[],'Spearmanr mean[95CI]':[],'AUROC mean[95CI]':[], 'AUPRC mean[95CI]':[]}
    
    pred_all = []
    gs_all = []
    
    for i, (train_idx,test_idx) in enumerate(kf.split(df_invivo)):
        print('Start preparing fold', i, '...')
        path = './invivo/fold_'+str(i)
        os.makedirs(path, exist_ok = True)
        TRAIN = df_invivo.iloc[train_idx]
        TEST = df_invivo.iloc[test_idx]
        TRAIN.to_csv(path+"/Train.csv", index = False)
        TEST.to_csv(path+"/Test.csv", index = False)
        TRAIN_data = TRAIN.iloc[:,4:]
        TEST_data = TEST.iloc[:,4:]
        
        # train in vivo model
        if model_type == 'lgb':
            predictor = train_lighgbm_model(TRAIN_data)
        elif model_type == 'xgb':
            predictor = train_xgboost_model(TRAIN_data)
        elif model_type == 'rf':
            predictor = train_rf_model(TRAIN_data)
        elif model_type == 'lr':
            predictor = train_lr_model(TRAIN_data)
        elif model_type == 'gpr':
            predictor = train_gpr_model(TRAIN_data)
        else:
            print('No model type:', model_type)

        filename = './params/fold_'+str(i)+'_model.sav'
        print('Saving lighgbm model ...')
        pickle.dump(predictor, open(filename, 'wb'))
        
        #predict on test set
        est=pickle.load(open(filename, 'rb'))
        pred=est.predict(TEST_data.iloc[:,:-1])
        pred = np.array(list(pred))
        gs = np.array(list(TEST_data.iloc[:,-1]))
        
        pred_all.extend(list(pred))
        gs_all.extend(list(gs))

        auroc = compute_auroc(pred,gs)
        auprc = compute_auprc(pred,gs)
        spearman_cor, _ = spearmanr(pred,gs)
        pearson_cor, _ = pearsonr(pred,gs)
        print("AUROC =", auroc, "AUPRC =", auprc,"Pearson's r =",pearson_cor, "Spearman's r =", spearman_cor)
        # save evaluations
        eva_df['fold'].append(str(i))
        eva_df['AUROC'].append(auroc)
        eva_df['AUPRC'].append(auprc)
        eva_df['Pearsonr'].append(pearson_cor)
        eva_df['Spearmanr'].append(spearman_cor)

    # Overall confidence analysis from k-fold results
    ci = 0.95
            
    mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, pearsonr_cor, ci)
    print("Mean[%d%sCI] Pearson's correlation is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
    eva_conf_df['Pearsonr mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))
            
    mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, spearmanr_cor, ci)
    print("Mean[%d%sCI] Spearman's correlation is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
    eva_conf_df['Spearmanr mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))
            
    mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, compute_auroc, ci)
    print("Mean[%d%sCI] AUROC is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
    eva_conf_df['AUROC mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))

    mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, compute_auprc, ci)
    print("Mean[%d%sCI] AUPRC is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
    eva_conf_df['AUPRC mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))

    eva_df = pd.DataFrame.from_dict(eva_df)
    eva_conf_df = pd.DataFrame.from_dict(eva_conf_df)

    return eva_df, eva_conf_df

def transfer_test_on_in_vitro(df_invitro):
    """ Transfer validation on in vitro dataset
    
    Parameters:
    -----------
    df_in_vitro: a Pandas dataframe
        preprocessed in vitro dataset
    
    Yields:
    -------
    eva_df: a Pandas dataframe
        performance of five in vivo models on in vitro dataset
    eva_conf_df: a Pandas dataframe
        confidence evaluations of in vivo models on in vitro dataset
    """
    os.makedirs('./invitro/', exist_ok = True)
    # split transfer testing set by 'Timepoint'(6HR/24HR) and 'Treatment'(DHA/UT)
    eva_df = {'data':[],'fold':[],'Pearsonr':[],'Spearmanr':[],'C-index':[]}
    eva_conf_df = {'data':[],'Pearsonr mean[95CI]':[],'Spearmanr mean[95CI]':[],'C-index mean[95CI]':[]}
    for i in sorted(set(df_invitro['Timepoint'])):
        for j in sorted(set(df_invitro['Treatment'])):
            path = './invitro/'+i+'_'+j
            data = i+'_'+j
            print(data)
            os.makedirs(path, exist_ok = True)
            TRANSFER = df_invitro.loc[(df_invitro['Timepoint'] == i)&(df_invitro['Treatment'] == j)] #select corresponding datasets
            TRANSFER_data = TRANSFER.iloc[:,4:]
            TRANSFER.to_csv(path+"/Test.csv", index = False)
            pred_all = []
            gs_all = []
            for k in range(10):
                filename = './params/fold_'+str(k)+'_model.sav'
                est=pickle.load(open(filename, 'rb'))
                predictions=est.predict(TRANSFER_data.iloc[:,:-1])
                TRANSFER['pred'] = predictions
                # compute mean prediction for each isolate
                pred = []
                gs = []
                for isolate in sorted(list(set(TRANSFER['Isolate']))):
                    pred.append(np.mean(TRANSFER.loc[TRANSFER['Isolate'] == isolate,'pred'].to_list()))
                    gs.append(np.mean(TRANSFER.loc[TRANSFER['Isolate'] == isolate,'DHA_IC50'].to_list()))
                
                pred_all.extend(pred)
                gs_all.extend(gs)
                # evaluations of predictions
                spearman_cor, _ = spearmanr(pred,gs)
                pearson_cor, _ = pearsonr(pred,gs)
                cidx = c_index(pred,gs)
                print("Pearson's r",pearson_cor,"Spearman's r = ", spearman_cor, "C-idex =", cidx)
                eva_df['data'].append(data)
                eva_df['fold'].append(k)
                eva_df['Pearsonr'].append(pearson_cor)
                eva_df['Spearmanr'].append(spearman_cor)
                eva_df['C-index'].append(cidx)
    
            # Overall confidence analysis from k-fold results
            ci = 0.95
            
            eva_conf_df['data'].append(data)
            
            mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, pearsonr_cor, ci)
            print("Mean[%d%sCI] Pearson's correlation is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
            eva_conf_df['Pearsonr mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))
            
            mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, spearmanr_cor, ci)
            print("Mean[%d%sCI] Spearman's correlation is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
            eva_conf_df['Spearmanr mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))
            
            mb, lb, ub = boostrapping_confidence_interval(pred_all, gs_all, c_index, ci)
            print("Mean[%d%sCI] C-idex is: %.4f[%.4f, %.4f]" % (int(ci*100), '%', mb, lb, ub))
            eva_conf_df['C-index mean[95CI]'].append("%.4f[%.4f, %.4f]" %(mb, lb, ub))
    
    eva_df = pd.DataFrame.from_dict(eva_df)
    eva_conf_df = pd.DataFrame.from_dict(eva_conf_df)
    return eva_df, eva_conf_df
