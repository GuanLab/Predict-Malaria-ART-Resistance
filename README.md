# Predict-Malaria-ART-Resistance

Winning algorithm in [2019 Malaria DREAM Challenge SubChallenge 2](https://www.synapse.org/#!Synapse:syn16924919/wiki/583955).

Transfer learning from in vivo to in vitro samples.

## Usage

```
python main.py [-h] [--train_path TRAIN_PATH] [--valid_path VALID_PATH]
               [-m MODEL_TYPE] [--no_quantile] [-n USE_TOP_FEATURES]
               [--feature_selection_mode]

Pipeline for buidling Malaria Artemisinin resistance in vivo prediction models and transfer learning on in vitro datasets.

optional arguments:
  -h, --help            show this help message and exit
  --train_path TRAIN_PATH
                        path to your training data, in .csv format
  --valid_path VALID_PATH
                        path to your transfer validation data, in .csv format
  -m MODEL_TYPE, --model_type MODEL_TYPE
                        machine learning models to use:
                                    lgb: LightGBM; 
                                    xgb: XGBoost;
                                    rf: Random Forest; 
                                    gpr: Gaussian Process Regression;
                                    lr: Linear Regression;
                                    default: lgb
  --no_quantile         if specified, do not use quantile normalization.
  -n USE_TOP_FEATURES, --use_top_features USE_TOP_FEATURES
                        if specified, used top features based on shap analysis.
                                    for example: 3, 5, 10, 20 or 30; 
                                    Otherwise, use the whole genome
  --feature_selection_mode
                        if used, begins leave-one-out feature selection strategy from specified gene set

```

It will generate the following folders:

`./invivo/`: preprocessed in vivo datasets  for model training and k-fold cross validation

`./invitro/`: in vitro dataset for transfering test

`./params/`: trained machine learning model parameters

`./performance/`: model performance in k-fold cross validation and transfering test

