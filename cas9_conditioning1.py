# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:34:27 2017

Analysis of conditioning data (first round) for cas9

@author: LocalAdmin1
"""

# Need to run cas9_assemble data first



def condhistFig(data, msk, color='grey', meancolor='black'):
    plt.figure()
    a = data[msk]
    for i in a:
        plt.plot(i, c=color, alpha=0.4)
    
    mean = np.mean(a)   
    plt.plot(mean, c=meancolor, linewidth=2)
    
    return mean

def cond2Dfig(df1, df2, key, sol='casein', ax=ax):
    dietmsk = df1.diet == 'np'
   
    a = [[df1[key][dietmsk], df2[key][dietmsk]],
          [df1[key][~dietmsk], df2[key][~dietmsk]]]

    x = data2obj2D(a)
    
    if sol == 'casein':
        barfacecolor = [col['np_cas'], col['np_cas'], col['lp_cas'], col['lp_cas']]
    else:
        barfacecolor = [col['np_malt'], col['np_malt'], col['lp_malt'], col['lp_malt']]
        
    ax, x, _, _ = jmfig.barscatter(x, paired=True,
                 barfacecoloroption = 'individual',
                 barfacecolor = barfacecolor,
                 scatteredgecolor = ['xkcd:charcoal'],
                 scatterlinecolor = 'xkcd:charcoal',
                 grouplabel=['NR', 'PR'],
                 scattersize = 60,
                 ax=ax)

# Analysis of conditioning Data
    
dfc1 = pd.DataFrame([(rats[x].casein1) for x in rats])
dfc2 = pd.DataFrame([(rats[x].casein2) for x in rats])
dfm1 = pd.DataFrame([(rats[x].maltodextrin1) for x in rats])
dfm2 = pd.DataFrame([(rats[x].maltodextrin2) for x in rats])

for df in [dfc1, dfc2, dfm1, dfm2]:
    df.insert(0,'ratid', [x for x in rats])
    df.insert(1,'diet', [rats[x].diet for x in rats])

df = pd.concat([dfc1, dfc2, dfm1, dfm2])

df.insert(2,'sol',['c']*48 + ['m']*48)
df.insert(3,'day',['c1']*24 + ['c2']*24 + ['m1']*24 + ['m2']*24)

df2 = df[['ratid', 'diet']][:48]
casall = df[df['day'] == 'c1']['total'] + df[df['day'] == 'c2']['total']
maltall = df[df['day'] == 'm1']['total'] + df[df['day'] == 'm2']['total']

df2.insert(2,'sol',['c']*24 + ['m']*24)
df2.insert(3,'total', casall.append(maltall))

# Fig showing total consumption across conditioning days               
fig = plt.figure(figsize=(3.2, 2.4))
ax = plt.subplot(1,1,1)
nplp2Dfig(df2, 'total', ax=ax)
plt.yticks([0, 5000, 10000, 15000], ('0', '5', '10', '15'))
ax.set_ylabel('Licks (x1000)')
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\03_condlicks.eps')
plt.title('Consumption during conditioning')


# Fig showing licking with each conditioning daya separately
fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(6.4, 2.4))
cond2Dfig(dfc1, dfc2, 'total', sol='casein', ax=ax[0])
cond2Dfig(dfm1, dfm2, 'total', sol='maltodextrin', ax=ax[1])
plt.title('Consumption during conditioning over both days')

# Fig showing palatability across conditioning days
fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(6.4, 2.4))
cond2Dfig(dfc1, dfc2, 'bMean', sol='casein', ax=ax[0])
cond2Dfig(dfm1, dfm2, 'bMean', sol='maltodextrin', ax=ax[1])
plt.title('Average cluster size during conditioning over both days')

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(6.4, 2.4))
cond2Dfig(dfc1, dfc2, 'bNum', sol='casein', ax=ax[0])
cond2Dfig(dfm1, dfm2, 'bNum', sol='maltodextrin', ax=ax[1])
plt.title('Average cluster number during conditioning over both days')

# Stats with R
r_df = df2[['ratid', 'sol', 'diet', 'total']]
ro.globalenv['r_df'] = r_df

ro.r('condlicks = aov(formula = total ~ sol * diet + Error(ratid / sol), data = r_df)')

print('Licks during conditioning')
print(ro.r('summary(condlicks)'))


## Fig showing histograms of consumption across conditioning days
#
#np_c1_x = condhistFig(data_cond['hist_cas1'], msk)
#np_c2_x = condhistFig(data_cond['hist_cas2'], msk)
#np_m1_x = condhistFig(data_cond['hist_malt1'], msk)
#np_m2_x = condhistFig(data_cond['hist_malt2'], msk)
#
#lp_c1_x = condhistFig(data_cond['hist_cas1'], ~msk)
#lp_c2_x = condhistFig(data_cond['hist_cas2'], ~msk)
#lp_m1_x = condhistFig(data_cond['hist_malt1'], ~msk)
#lp_m2_x = condhistFig(data_cond['hist_malt2'], ~msk)
