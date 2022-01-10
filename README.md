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

## External validation datasets besides DREAM challenges

While our model was designed to be applied on a challenging task to predict ART resistance based on transcripton profiles sample at different conditions (in vivo, in vitro, ex vivo) and different microarray platforms in DREAM challenge, we also tested our model on other independent public datasets for predicting Malaria ART resistance, which is put under `./exteral_data`.

* in vivo: GSE149735 (Zhu et al. 2021) *ps: this data is not included since correct label is not yet available from [paper](https://www.biorxiv.org/content/10.1101/2021.05.17.444396v1)*
* in vitro: GSE151189 GSE61536 
* ex vivo: GSE25878 GSE59098


first install GEOparse:
`pip install GEOparse`

then pull down public transcriptomes:
```
cd ./exteral_data
python getGEO.py
```
This will both pull down the soft files from GEO and create `.tsv` files for test.

```
ex_vitro_GSE25878.tsv  
ex_vitro_GSE59098.tsv
in_vitro_GSE151189.tsv  
in_vitro_GSE61536.tsv
```




