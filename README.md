# Predict-Malaria-ART-Resistance

Winning algorithem in [2019 Malaria DREAM Challenge SubChallenge 2](https://www.synapse.org/#!Synapse:syn16924919/wiki/583955).

## Transfer learning between in vivo and in vitro samples

By running
```
python main.py <model_type> 
```
it will generate the following folders:

`./invivo/`: preprocessed in vivo datasets  for model training and k-fold cross validation

`./invitro/`: in vitro dataset for transfering test

`./params/`: trained machine learning model parameters

`./performance/`: model performance in k-fold cross validation and transfering test

