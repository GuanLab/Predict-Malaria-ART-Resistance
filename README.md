# Predict-Malaria-ART-Resistance
## Train Malaria ART resistance prediction model
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

