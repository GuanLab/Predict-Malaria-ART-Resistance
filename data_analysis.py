import pandas as pd
from collections import Counter
for i in range(10):
    p = pd.read_csv('./invivo/fold_'+str(i)+'/Train.csv')
    print('fold_',i)
    print(dict(Counter(p['Kmeans.Grp'])))
    print(dict(Counter(p['Country'])))
    print(dict(Counter(p['Asexual.stage..hpi.'])))
    #print(dict(Counter(p['ClearanceRate'])))
