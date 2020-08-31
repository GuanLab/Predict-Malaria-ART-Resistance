import numpy as np
import sklearn.ensemble
# light GBM modules
import lightgbm as lgb
# XGboost modules
from xgboost import XGBRegressor
# GPR modules
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
# linear regression modules
from sklearn.linear_model import LinearRegression
# random forest modules
import sklearn.ensemble

""" All types of model tested in this study
"""

def train_lighgbm_model(TRAIN):
    """ Train LightGBM Regression model
    Params:
    -------
    TRAIN: Pandas dataframe
        training set; last column is label
    Yields:
    -------
    lgb: LightGBM model
    """
    train_dat=TRAIN.iloc[:,:-1]
    train_gs=TRAIN.iloc[:,-1]
    lgb_train = lgb.Dataset(np.asarray(train_dat), np.asarray(train_gs))
    params = {
            'boosting_type': 'gbdt',
            'objective': 'regression',
            'num_leaves': 5,
            'learning_rate': 0.05,
            'verbose': 0,
            'n_estimators': 800,
            'reg_alpha': 2.0,
            #'first_metric_only': True
            }
    gbm = lgb.train(params,
            train_set = lgb_train,
            #valid_sets = lgb_test,
            #early_stopping_rounds = 5000,
            num_boost_round=1000
            )
    return gbm

def train_xgboost_model(TRAIN):
    """ Train XGboost Regression model
    Params:
    -------
    TRAIN: Pandas dataframe
        training set; last column is label
    Yields:
    -------
    xgp: XGboost Regression model
    """
    train_dat=TRAIN.iloc[:,:-1]
    train_gs=TRAIN.iloc[:,-1]
    xgb = XGBRegressor()
    xgb.fit(train_dat, train_gs)
    return xgb

def train_rf_model(TRAIN):
    """ Train Random Forest Regression model
    Params:
    -------
    TRAIN: Pandas dataframe
        training set; last column is label
    Yields:
    -------
    clf: Random Forest Regression model
    """
    train_dat=TRAIN.iloc[:,:-1]
    train_gs=TRAIN.iloc[:,-1]
    clf = sklearn.ensemble.RandomForestRegressor(n_estimators=100, max_depth=2, random_state=0)
    clf.fit(train_dat, train_gs)
    return clf


def train_gpr_model(TRAIN):
    """ Train Gaussian Process Regression model
    Params:
    -------
    TRAIN: Pandas dataframe
        training set; last column is label
    
    Yields:
    -------
    gpr: Gaussian Process Regression model
    """
    train_dat=TRAIN.iloc[:,:-1]
    train_gs=TRAIN.iloc[:,-1]
    kernel = DotProduct() + WhiteKernel()
    gpr = GaussianProcessRegressor(kernel=kernel, random_state=0).fit(train_dat, train_gs)
    return gpr

def train_lr_model(TRAIN):
    """ Train Linear Regression model
    Params:
    -------
    TRAIN: Pandas dataframe
        training set; last column is label
    
    Yields:
    -------
    lr: linear regression model
    """
    train_dat=TRAIN.iloc[:,:-1]
    train_gs=TRAIN.iloc[:,-1]
    lr =  LinearRegression().fit(train_dat, train_gs)
    return lr
