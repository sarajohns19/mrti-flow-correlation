import os
import pandas as pd
import numpy as np
import openpyxl as xl

# import .csv as dataframe
#datapath = '/Users/sarajohnson/Box Sync/Flow Cytometry Data/Compiled Results/'
datapath = '/Users/sarajohnson/Box Sync/Work Documents/XFUS Research/Immunotherapy Project/Experiments/'
wb = xl.load_workbook(filename= datapath + 'Growth_Data_Import.xlsx')

ws = wb['Days Since Treatment']
subject = []
days_since =[]
rowlist = list(ws.values)
for row in rowlist[1:]:
    subject.append(row[0])
    days_since.append([x for x in row[1:] if x is not None])

ws = wb['Volume']
volume =[]
rowlist = list(ws.values)
for row in rowlist[1:]:
    volume.append([x for x in row[1:] if x is not None])

ws = wb['% Change in Volume']
volume_change =[]
rowlist = list(ws.values)
for row in rowlist[1:]:
   volume_change.append([x for x in row[1:] if x is not None])


# create dataframe

tumor_volume_df = pd.DataFrame.from_dict({'subject': subject, 'Days of Treatment': days_since,
                                         'Tumor Volume': volume, 'Tumor Volume Change': volume_change})

#%%

tumor_volume_df['Start Volume'] = [x[0] for x in tumor_volume_df['Tumor Volume']]
tumor_volume_df['End Volume'] = [x[-1] for x in tumor_volume_df['Tumor Volume']]
tumor_volume_df['End % Change in Volume'] = [x[-1] for x in tumor_volume_df['Tumor Volume Change']]

#%% save dataframe

tumor_volume_df.to_pickle('tumor_volume_data.pkl')