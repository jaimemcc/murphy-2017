# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 14:48:26 2019

@author: admin
"""

# need to run cas9_assembledata.py first

import pandas as pd

pref_sessions = {}
for session in sessions:
    x = sessions[session]
    if x.sessioncode == 'preference1':
        pref_sessions[x.sessionID] = x


df = pd.DataFrame([sessions[s].rat for s in pref_sessions], columns=['rat'])
df['session'] = [sessions[s].session for s in pref_sessions]

df['diet'] = [sessions[s].diet for s in pref_sessions]

df['cas_alpha'] = [sessions[s].lickdata_cas['weib_alpha'] for s in pref_sessions]
df['cas_beta'] = [sessions[s].lickdata_cas['weib_beta'] for s in pref_sessions]
df['cas_rsq'] = [sessions[s].lickdata_cas['weib_rsq'] for s in pref_sessions]

df['malt_alpha'] = [sessions[s].lickdata_malt['weib_alpha'] for s in pref_sessions]
df['malt_beta'] = [sessions[s].lickdata_malt['weib_beta'] for s in pref_sessions]
df['malt_rsq'] = [sessions[s].lickdata_malt['weib_rsq'] for s in pref_sessions]


def nplp2Dfig(df, key_cas, key_malt, ax):
    dietmsk = df.diet == 'np'
    
    a = [[df[key_cas][dietmsk], df[key_malt][dietmsk]],
          [df[key_cas][~dietmsk], df[key_malt][~dietmsk]]]

    x = data2obj2D(a)

    ax, x, _, _ = jmfig.barscatter(x, paired=True,
                 barfacecoloroption = 'individual',
                 barfacecolor = [col['np_cas'], col['np_malt'], col['lp_cas'], col['lp_malt']],
                 scatteredgecolor = ['xkcd:charcoal'],
                 scatterlinecolor = 'xkcd:charcoal',
                 grouplabel=['NR', 'PR'],
                 scattersize = 60,
                 ax=ax)

f, ax = plt.subplots(figsize=(8, 3), ncols=3)
nplp2Dfig(df, 'cas_alpha', 'malt_alpha', ax[0])

nplp2Dfig(df, 'cas_beta', 'malt_beta', ax[1])

nplp2Dfig(df, 'cas_rsq', 'malt_rsq', ax[2])