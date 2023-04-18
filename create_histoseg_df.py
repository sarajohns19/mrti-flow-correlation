import os
import pandas as pd
import numpy as np

datapath = '/Users/sarajohnson/Box Sync/Work Documents/XFUS Research/Immunotherapy Project/Experiments/'
histo_df = pd.read_csv(datapath + 'H&E_segmentations.csv')

histo_df.rename(columns=lambda x: x.split(': ')[-1], inplace=True)
histo_df.rename(columns={'Image': 'subject', 'Necrosis %': 'HE_Necrosis1', 'Negative %': 'HE_Negative',
                        'Pinky %': 'HE_Necrosis2', 'live cell %': 'HE_Live'}, inplace=True)


histo_df['subject']= [int(x[x.find('_M')+2:x.find('.')]) for x in histo_df['subject']]
histo_df[histo_df.columns[1:]] = histo_df[histo_df.columns[1:]]/100
histo_df.to_pickle('h&e_necrosis_seg.pkl')

