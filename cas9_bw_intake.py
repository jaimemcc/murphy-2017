# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:33:40 2017

Analysis of food and body weight data for cas9

@author: LocalAdmin1
"""

import numpy as np
import os

import pandas as pd
import matplotlib.pyplot as plt

import rpy2.robjects as ro

from rpy2.robjects import r, pandas2ri, numpy2ri
pandas2ri.activate()
numpy2ri.activate()

from scipy import stats

import JM_custom_figs as jmfig

userhome = os.path.expanduser('~')
metafile = userhome + '\\Dropbox\\Python\\cas9\\CAS9bw_metafile.csv'

data = pd.read_csv(metafile, index_col=['rat', 'diet'])

# Statistics

data = data[:].stack()
data = data.to_frame()
data.reset_index(inplace=True) 
data.columns = ['rat', 'diet', 'day', 'licks']
ro.globalenv['r_data'] = data

ro.r('bodyweight = aov(formula = licks ~ day * diet + Error(rat / day), data = r_data)')

print(ro.r('summary(bodyweight)'))


data = pd.read_csv(metafile, index_col=['rat'])

np_mean = data[data['diet'] == 'np'].mean()
np_sem = data[data['diet'] == 'np'].std() / np.sqrt(len(data['diet'] == 'np'))

lp_mean = data[data['diet'] == 'lp'].mean()
lp_sem = data[data['diet'] == 'lp'].std() / np.sqrt(len(data['diet'] == 'lp'))

# Figure 1A - Body weight
fig = plt.figure(figsize=(3.2,2.4))
ax = plt.subplot(1,1,1)
np_mean.plot(yerr=np_sem, color='xkcd:charcoal', marker='o', markerfacecolor='white')
lp_mean.plot(yerr=lp_sem, color='xkcd:kelly green', marker='o', markerfacecolor='white')
ax.set_ylim([400, 550])
ax.set_xlim([-1, 17])
plt.xticks([1,6,11,16], ('0', '5', '10', '15'))
plt.yticks([400, 450, 500, 550])
ax.set_ylabel('Body weight (g)')
ax.set_xlabel('Days since diet switch')

plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\01_bodyweight.eps')


# Figure 1B - Food intake
# Data from CAS9_exptdetails.xls > fi_cas9+cas56

foodintake_np = [22.0, 25.6, 27.6, 26.1]
foodintake_lp = [25.5, 27.3, 29.1, 28.5]

foodintake_npAll = [22.0, 25.6, 27.6, 26.1, 22.9, 21.0, 21.3, 20.6]
foodintake_lpAll = [25.5, 27.3, 29.1, 28.5, 27.4, 29.9, 23.5, 26.3]

fi = data2obj1D([foodintake_np, foodintake_lp])
#fi = np.array([foodintake_np, foodintake_lp],
#         dtype='object')
mpl.rcParams['figure.subplot.left'] = 0.25
fig = plt.figure(figsize=(1.8,2.4))
ax = plt.subplot(1,1,1)
jmfig.barscatter(fi, barfacecoloroption='individual',
                 barwidth = 0.8,
                 barfacecolor = ['xkcd:silver', 'xkcd:kelly green'],
                 scatteredgecolor = ['xkcd:charcoal'],
                 scattersize = 40,
                 ylabel = 'Average food intake (g/day)',
                 grouplabel=['NR', 'PR'],
                 ax=ax)
plt.yticks([0, 10, 20, 30])
ax.set_xlim([0.25,2.75])
ax.set_ylim([0, 35])
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\02_foodintake.eps')

fi = data2obj1D([foodintake_npAll, foodintake_lpAll])
fig = plt.figure(figsize=(1.8,2.4))
ax = plt.subplot(1,1,1)
jmfig.barscatter(fi, barfacecoloroption='individual',
                 barwidth = 0.8,
                 barfacecolor = ['xkcd:silver', 'xkcd:kelly green'],
                 scatteredgecolor = ['xkcd:charcoal'],
                 scattersize = 40,
                 grouplabel=['NR', 'PR'],
                 ax=ax)
plt.yticks([0, 10, 20, 30])
ax.set_xlim([0.25,2.75])
ax.set_ylim([0, 35])
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\02b_foodintake.eps')

mpl.rcParams['figure.subplot.left'] = 0.15
fi_stats = stats.ttest_ind(foodintake_npAll, foodintake_lpAll)
print(fi_stats)