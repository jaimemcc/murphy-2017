# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 11:39:15 2017

This script analyzes conditioning and preference data together.
Requires cas9_assemble data, conditioning1, and pref1 to be run first

@author: James Rig
"""

def condprefFig(cdata, mdata, key):
    
    fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(6.4, 2.4))
    
    x = [1,2,3]
    
    colors=[]
    colors.append(['xkcd:silver', 'white'])
    colors.append(['xkcd:kelly green', 'xkcd:light green'])
    
    for diet, axis in zip(['np', 'lp'], [0, 1]):
        msk = cdata[0].diet == diet
        casMean = [x[key][msk].mean() for x in cdata]
        casSEM = [x[key][msk].std()/np.sqrt(sum(msk)) for x in cdata]
        
        maltMean = [x[key][msk].mean() for x in mdata]
        maltSEM = [x[key][msk].std()/np.sqrt(sum(msk)) for x in mdata]
        
        for avg, sem, cols in zip([casMean, maltMean], [casSEM, maltSEM], [0, 1]):
            ax[axis].errorbar(x, avg, yerr=sem,
                              marker='o',
                              markersize='8',
                              color='xkcd:charcoal',
                              markerfacecolor=colors[axis][cols],
                              markeredgecolor='xkcd:charcoal')
#        ax[axis].errorbar(x, maltMean, yerr=maltSEM, marker='o', mec=colors[axis][1])
        
    ax[axis].set_xlim([0.5, 3.5])
    ax[0].set_ylim([0, 60])
    plt.xticks([1,2,3],['Cond. 1','Cond. 2','Pref.'])

    return fig, ax

fig, ax = condprefFig([dfc1, dfc2, dfc], [dfm1, dfm2, dfm], 'bMean')
ax[0].set_ylabel('Mean bout size')

plt.title('Palatability across conditioning and preference test')



# Analysis of palatability over conditioning and test
#key = 'bMean'
#msk = dfc1.diet == 'lp'
#
#z = dfc1[key][msk].mean()
#z2 = dfc2[key][msk].mean()
#z3 = dfc[key][msk].mean()
#
#y = dfm1[key][msk].mean()
#y2 = dfm2[key][msk].mean()
#y3 = dfm[key][msk].mean()


