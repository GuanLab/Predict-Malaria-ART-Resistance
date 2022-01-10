import GEOparse
import pandas as pd
# in vivo: GSE59099(DREAM training set, Mok et al. 2015); GSE149735 (Zhu et al. 2021) !!correct label is not available from paper
# in vitro: GSE151189 GSE62132;GSE61536 
# ex vivo: GSE25878 GSE59098


# in vivo
gse = GEOparse.get_GEO(geo="GSE149735", destdir="./")

gpl_name, gpl = list(gse.gpls.items())[0]
plt_orf = gpl.table[['ORF', 'ORF_old']].drop_duplicates(subset = None).reset_index(drop = True)
# map old orf id to new id
orf_map = {i['ORF_old']:i['ORF'] for _,i in plt_orf.iterrows()}

"""
    # GSE149735
    # treatment condition
    t = gsm.metadata['title'] #T20002_0hour_0day
    acc = gsm.metadata['geo_accession'] #GSM4510099
    pl = gsm.metadata['platform_id'] #GPL18893
    des = gsm.metadata['description'] #KH001-2-010_H00D00 # is this patient id?
    # map the patient id to clearnance half life
    ch1 = gsm.metadata['characteristics_ch1'] #['sample type/stage: In vivo, Blood stage', 'treatment: Prior to ACT treatment, 0h', 'batch: 1']
    ch1[0]  #sample type/stage: In vivo, Blood stage'
    ch1[1]  #'treatment: Prior to ACT treatment, 0h'
    ch1[2]  #'batch: 17'
"""



"""

# in vitro
gse = GEOparse.get_GEO(geo="GSE151189", destdir="./") # PMID: 33483501

# platform gse.gpls 
assert len(gse.gpls) == 1, "Platform more than 1!"

# platform name and platform table
gpl_name, gpl = list(gse.gpls.items())[0]

print("Platform name:", gpl_name)
print("Platform table:", gpl.table)

# platform in pandas dataframe format
plt = gpl.table

# samples gse.gsms

all_title = []
all_exp = []
all_label = []  # resistance, 1 for  resistant, 0 for non-resistant
for gsm_name, gsm in gse.gsms.items():
    '''
    # GSE151189
    '''
    t = gsm.metadata['title'][0] # Dd2_[WT|R539T|C580Y|]_[0|3|6|16|24|48]h_rep[1|2|3] , Cam3II_[revWT|R539T|]_[0|8|16|24|32|40|48]h_rep[1|2|3]
    print(t)
    '''
    strain:
        Dd2: lab strain
        Cam3II: isolate from patients
    mutation:
        R539T: resistant ++
        C580Y: resistant +
        WT: resistant -
        *for resistance we denote both mutant as 1 and WT as 0
    0-48h:
        hours after DHA treatment
    rep1-3:
        biological replicates
    '''
    exp =gsm.table.rename(columns = {'ID_REF':'ID'})
    # map probe to gene by left join
    new_exp = exp.merge(plt, on='ID', how =  'left')[['ID', 'ORF', 'VALUE']].groupby('ORF').mean()[['VALUE']].T
    if t.split('_')[1] in ['R539T', 'C580Y']:
        all_label.append(1)
    else:
        all_label.append(0)
    all_title.append(t)
    all_exp.append(new_exp)

df_exp= pd.concat(all_exp).reset_index(drop = True)
df_exp['sample'] = all_title
df_exp['label'] = all_label
df_exp.to_csv('in_vitro_GSE151189.csv')

"""

#in vitro
gse = GEOparse.get_GEO(geo="GSE61536", destdir="./") # PMID: 26490244

# platform gse.gpls 
assert len(gse.gpls) == 1, "Platform more than 1!"

# platform name and platform table
gpl_name, gpl = list(gse.gpls.items())[0]

print("Platform name:", gpl_name)
print("Platform table:", gpl.table)

# platform in pandas dataframe format
plt = gpl.table

# samples gse.gsms

all_title = []
all_exp = []
all_label = []  # resistance, 1 for  resistant, 0 for non-resistant

for gsm_name, gsm in gse.gsms.items():
    t = gsm.metadata['title'][0] # [dha|vehicle control] 1 h rep 1 of 5
    print(t)
    ''' 
    strain: Pf K1 strain, a lab ART resistant strain
    deve_stage: ring|troph|schiz
    we only use DHA treatment in this study for testing 
        hours after DHA treatment
    rep1-5:
        biological replicates
    '''
    if t.startswith('dha'):
        print(t)
        exp =gsm.table.rename(columns = {'ID_REF':'ID'})
        # map probe to gene by left join
        new_exp = exp.merge(plt, on='ID', how =  'left')[['ID', 'ORF', 'VALUE']].groupby('ORF').mean()[['VALUE']].T.rename(columns= orf_map)
        all_title.append(t)
        all_exp.append(new_exp)
        all_label.append(1)

df_exp= pd.concat(all_exp).reset_index(drop = True)
df_exp['sample'] = all_title
df_exp['label'] = all_label
df_exp.to_csv('in_vitro_GSE61536.csv')

"""
# NOT avaible since it's RNA seq, not microarray
gse = GEOparse.get_GEO(geo="GSE62132", destdir="./") # PMID: 26490244

# platform gse.gpls
assert len(gse.gpls) == 1, "Platform more than 1!"

# platform name and platform table
gpl_name, gpl = list(gse.gpls.items())[0]

print("Platform name:", gpl_name)
print("Platform table:", gpl.table)

# platform in pandas dataframe format
plt = gpl.table

# samples gse.gsms

all_title = []
all_exp = []
all_label = []  # resistance, 1 for  resistant, 0 for non-resistant

for gsm_name, gsm in gse.gsms.items():
    t = gsm.metadata['title'][0] # Pf_[ring|troph|schiz]_[un|vehicle|DHA]_rep[1|2]
    print(t)
    '''
    strain: Pf K1 strain, a lab ART resistant strain
    deve_stage: ring|troph|schiz
    we only use DHA treatment in this study for testing
        hours after DHA treatment
    rep1-2:
        biological replicates
    '''
    if t.split('_')[2] =='DHA':
        exp =gsm.table.rename(columns = {'ID_REF':'ID'})
        # map probe to gene by left join
        new_exp = exp.merge(plt, on='ID', how =  'left')[['ID', 'ORF', 'VALUE']].groupby('ORF').mean().T
        all_title.append(t)
        all_exp.append(new_exp)
        all_label.eppend(1)

df_exp= pd.concat(all_exp).reset_index(drop = True)
df_exp['sample'] = all_title
df_exp['label'] = all_label
df_exp.to_csv('in_vitro_GSE62132.csv')
"""

# ex vivo
gse = GEOparse.get_GEO(geo="GSE25878", destdir="./") # PMID: 21810278

# platform gse.gpls
assert len(gse.gpls) == 1, "Platform more than 1!"

# platform name and platform table
gpl_name, gpl = list(gse.gpls.items())[0]

print("Platform name:", gpl_name)
print("Platform table:", gpl.table)

# platform in pandas dataframe format
plt = gpl.table

# samples gse.gsms

all_title = []
all_exp = []
all_label = []  # resistance, 1 for  resistant, 0 for non-resistant

for gsm_name, gsm in gse.gsms.items():
    t = gsm.metadata['title'][0] # ex-vivo_[( ART-sensitive)|( ART-resistant)]_[( Xepon Laos)|( Mae Sot Thailand)|( Pailin Cambodia)]_[BMT061]_[0-48]H
    print(t)
    '''
    all ex vivo datasets; 11 isolates collected from 3 demographic sites
    [0-48]H: the time duration of growing ex vivo 
    '''
    if t.split('_')[1]==' ART-resistance':
        all_label.append(1)
    else:
        all_label.append(0)

    exp =gsm.table.rename(columns = {'ID_REF':'ID'})
    # map probe to gene by left join
    new_exp = exp.merge(plt, on='ID', how =  'left')[['ID', 'ORF', 'VALUE']].groupby('ORF').mean()[['VALUE']].T.rename(columns= orf_map)
    all_title.append(t)
    all_exp.append(new_exp)

df_exp= pd.concat(all_exp).reset_index(drop = True)
df_exp['sample'] = all_title
df_exp['label'] = all_label
df_exp.to_csv('ex_vitro_GSE25878.csv')


# from the same isolates in Mok et al's paper in in vivo training set
# The samples with <=5 hours of clearance half life are labeled as  “Fast” in terms of clearance rate, and considered as non-ART-resistant samples
ref = pd.read_csv('../rawdata/mok_supplementary/Table3.Clearance_rate.csv')

gse = GEOparse.get_GEO(geo="GSE59098", destdir="./") # PMID: 21810278

# platform gse.gpls
assert len(gse.gpls) == 1, "Platform more than 1!"

# platform name and platform table
gpl_name, gpl = list(gse.gpls.items())[0]

print("Platform name:", gpl_name)
print("Platform table:", gpl.table)

# platform in pandas dataframe format
plt = gpl.table

# samples gse.gsms

all_title = []
all_exp = []
all_label = []  # resistance, 1 for  resistant, 0 for non-resistant

for gsm_name, gsm in gse.gsms.items():
    t = gsm.metadata['title'][0] # KH004-068-0h-312918
    print(t)
    '''
    all ex vivo datasets; 19 isolates collected from Pailin, Cambodia
    [0-48]h: the time duration of growing ex vivo 
    '''
    pid ='-'.join(t.split('-')[:2])
    try:
        clr = ref.loc[ref['Patient Code']==pid, 'parasite clearance halflife (h)'].to_list()[0]
        print(clr)
        if clr >=5:
            all_label.append(1)
        else:
            all_label.append(0)
        exp =gsm.table.rename(columns = {'ID_REF':'ID'})
        # map probe to gene by left join
        new_exp = exp.merge(plt, on='ID', how =  'left')[['ID', 'ORF', 'VALUE']].groupby('ORF').mean()[['VALUE']].T
        all_title.append(t)
        all_exp.append(new_exp)
    except:
        print(pid)

df_exp= pd.concat(all_exp).reset_index(drop = True)
df_exp['sample'] = all_title
df_exp['label'] = all_label
df_exp.to_csv('ex_vitro_GSE59098.csv')



