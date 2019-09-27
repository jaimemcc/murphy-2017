# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:37:32 2017

Analysis of preference day 1 for cas9

@author: LocalAdmin1
"""
# Need to run cas9_assembledata first

def prefhistFig(ax1, ax2, df, factor):
    casmsk = df.sol == 'c'
    dietmsk = df.diet == 'np'

    shadedError(ax1, df[factor][casmsk & dietmsk], linecolor='black')
    ax1 = shadedError(ax1, df[factor][~casmsk & dietmsk], linecolor='xkcd:bluish grey')
    ax1.set_xticks([0,10,20,30])
    ax1.set_xticklabels(['0', '20', '40', '60'])
    
    shadedError(ax2, df[factor][casmsk & ~dietmsk], linecolor='xkcd:kelly green')
    ax2 = shadedError(ax2, df[factor][~casmsk & ~dietmsk], linecolor='xkcd:light green')
    ax2.set_xticks([0,10,20,30])
    ax2.set_xticklabels(['0', '20', '40', '60'])
    
#    jmfig.setsameaxislimits([ax1, ax2])
    
def shadedError(ax, yarray, linecolor='black', errorcolor = 'xkcd:silver'):
    yarray = np.array(yarray)
    y = np.mean(yarray)
    yerror = np.std(yarray)/np.sqrt(len(yarray))
    x = np.arange(0, len(y))
    ax.plot(x, y, color=linecolor)
    ax.fill_between(x, y-yerror, y+yerror, color=errorcolor, alpha=0.4)
    
    return ax

def casVmaltFig(ax, df):
    # prepare data
    casmsk = df.sol == 'c'
    casdata = np.array(df['licks'][casmsk])
    maltdata = np.array(df['licks'][~casmsk])
    xydataAll = []
    for cas, malt in zip(casdata, maltdata):
        xydata = []
        x = np.array([cas[1:], [1]*(len(cas)-1)])
        y = np.array([malt[1:], [2]*(len(malt)-1)])
        alllicks = np.concatenate((x,y),axis=1)
        idx = np.argsort(alllicks[0])
        sortedlicks = alllicks[:,idx]
        xydata.append(np.cumsum(np.array(sortedlicks[1,:] == 1, dtype=int)))
        xydata.append(np.cumsum(np.array(sortedlicks[1,:] != 1, dtype=int)))
        xydataAll.append(xydata)
    
    dietmsk = (df.diet == 'np')
    dietmsk = dietmsk[:24]
    
    # plot line of unity    
    ax.plot([0, 5500], [0, 5500], '--', color='xkcd:silver', linewidth=0.5)
    
    npdata = [x for i,x in enumerate(xydataAll) if dietmsk[i]]
    for x in npdata:
        ax.plot(x[0], x[1], c='xkcd:silver', alpha=0.2, linewidth=0.5)
        ax.scatter(x[0][-1], x[1][-1], c='none', edgecolors='xkcd:charcoal')
    
    lpdata = [x for i,x in enumerate(xydataAll) if not dietmsk[i]]
    for x in lpdata:
        ax.plot(x[0], x[1], c='xkcd:light green', alpha=0.2, linewidth=0.5)
        ax.scatter(x[0][-1], x[1][-1], color='none', edgecolors='xkcd:kelly green')
        
    max_x = np.max([ax.get_xlim(), ax.get_ylim()])
    ax.set_xlim([-300, max_x])
    ax.set_ylim([-300, max_x])

    
    return xydataAll

# Analysis of Preference Day 1

dfc = pd.DataFrame([(rats[x].preference1_cas) for x in rats])
dfm = pd.DataFrame([(rats[x].preference1_malt) for x in rats])

for df in [dfc, dfm]:
    df.insert(0,'ratid', [x for x in rats])
    df.insert(1,'diet', [rats[x].diet for x in rats])
    
df = pd.concat([dfc, dfm])

df.insert(2,'sol',['c']*24 + ['m']*24)

df2 = df[['ratid', 'diet']][:24]
pref = df.total[:24]/(df.total[:24]+df.total[24:])
df2.insert(2,'pref', pref)

# Figure 3A - Licks over time, histogram
mpl.rcParams['figure.subplot.wspace'] = 0.1
mpl.rcParams['figure.subplot.left'] = 0.15
fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(3.2, 2.4))

prefhistFig(ax[0], ax[1], df, 'hist')
fig.text(0.55, 0.04, 'Time (min)', ha='center')
ax[0].set_ylabel('Licks per 2 min')
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\04_prefhist.eps')


# Figure 3B - casein licks vs maltodextrin licks
mpl.rcParams['figure.subplot.left'] = 0.25
fig = plt.figure(figsize=(2.4, 2.4))
ax = plt.subplot(1,1,1)                
casVmaltFig(ax, df)
ax.set_xlabel('Licks for casein')
ax.set_ylabel('Licks for maltodextrin')
plt.yticks([0, 2000, 4000, 6000])
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\05_casvmalt.eps')

mpl.rcParams['figure.subplot.left'] = 0.15

# Figure 3C - casein and malt preference (licks)

# Analysis of licks
fig = plt.figure(figsize=(3.2, 2.4))
ax = plt.subplot(1,1,1) 
nplp2Dfig(df, 'total', ax)
ax.set_ylabel('Licks (x1000)')
plt.yticks([0, 2000, 4000, 6000], ('0', '2', '4', '6'))
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\06_preftotal.eps')
plt.title('Casein vs. maltodextrin')

# Analysis of palatability

# Figure 4A - licks per burst

fig = plt.figure(figsize=(3.2, 2.4))
ax = plt.subplot(1,1,1)
nplp2Dfig(df, 'bMean', ax)
ax.set_ylabel('Average licks per cluster')
ax.set_yticks([0, 50, 100, 150])
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\07_prefburstmean.eps')
plt.title('Licks per burst')

fig = plt.figure(figsize=(3.2, 2.4))
ax = plt.subplot(1,1,1)
nplp2Dfig(df, 'bNum', ax)
ax.set_ylabel('Number of clusters')
ax.set_yticks([0, 50, 100, 150, 200])
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\08_prefburstnum.eps',
            transparent=True)
plt.title('Number of bursts')


# Figure 3D - protein preference

dietmsk = df2.diet == 'np'
a = data2obj1D([df2['pref'][dietmsk], df2['pref'][~dietmsk]])

mpl.rcParams['figure.subplot.left'] = 0.25
fig = plt.figure(figsize=(1.8, 2.4))
ax = plt.subplot(1,1,1)
jmfig.barscatter(a, barfacecoloroption = 'between', barfacecolor = ['xkcd:silver', 'xkcd:kelly green'],
                     scatteredgecolor = ['black'],
                     scatterlinecolor = 'black',
                     grouplabel=['NR', 'PR'],
                     barwidth = 0.8,
                     scattersize = 40,
                     ylabel = 'Casein preference',
                     ax=ax)
ax.set_yticks([0, 0.5, 1.0])
ax.set_xlim([0.25,2.75])
ax.set_ylim([0, 1.1])
ax.set_ylabel('Casein preference')
plt.savefig(userhome + '\\Dropbox\\Python\\cas9\\cas9_figs\\09_caseinpref.eps',
            transparent=True)
plt.title('Casein preference')

mpl.rcParams['figure.subplot.left'] = 0.15

## Statistics with R

r_df = df[['ratid', 'sol', 'diet', 'total', 'bMean', 'bNum']]
ro.globalenv['r_df'] = r_df

ro.r('totallicks = aov(formula = total ~ sol * diet + Error(ratid / sol), data = r_df)')
ro.r('burstMean = aov(formula = bMean ~ sol * diet + Error(ratid / sol), data = r_df)')
ro.r('burstNum = aov(formula = bNum ~ sol * diet + Error(ratid / sol), data = r_df)')

print('Total licks')
print(ro.r('summary(totallicks)'))

ro.r('np_casvmalt = t.test(r_df$total[r_df$diet=="np" & r_df$sol=="c"], r_df$total[r_df$diet=="np" & r_df$sol=="m"], paired=TRUE)')
print('Normal protein rats - casein vs. maltodextrin')
print(ro.r('np_casvmalt'))

ro.r('lp_casvmalt = t.test(r_df$total[r_df$diet=="lp" & r_df$sol=="c"], r_df$total[r_df$diet=="lp" & r_df$sol=="m"], paired=TRUE)')
print('LOW PROTEIN rats - casein vs. maltodextrin')
print(ro.r('lp_casvmalt'))

# Analysis of Licks per burst

print('Licks per burst')
print(ro.r('summary(burstMean)'))

ro.r('np_casvmalt = t.test(r_df$bMean[r_df$diet=="np" & r_df$sol=="c"], r_df$bMean[r_df$diet=="np" & r_df$sol=="m"], paired=TRUE)')
print('Normal protein rats (licks per burst) - casein vs. maltodextrin')
print(ro.r('np_casvmalt'))

ro.r('lp_casvmalt = t.test(r_df$bMean[r_df$diet=="lp" & r_df$sol=="c"], r_df$bMean[r_df$diet=="lp" & r_df$sol=="m"], paired=TRUE)')
print('LOW PROTEIN rats (licks per burst) - casein vs. maltodextrin')
print(ro.r('lp_casvmalt'))


print('Number of bursts')
print(ro.r('summary(burstNum)'))

ro.r('np_casvmalt = t.test(r_df$bNum[r_df$diet=="np" & r_df$sol=="c"], r_df$bNum[r_df$diet=="np" & r_df$sol=="m"], paired=TRUE)')
print('Normal protein rats (burst number) - casein vs. maltodextrin')
print(ro.r('np_casvmalt'))

ro.r('lp_casvmalt = t.test(r_df$bNum[r_df$diet=="lp" & r_df$sol=="c"], r_df$bNum[r_df$diet=="lp" & r_df$sol=="m"], paired=TRUE)')
print('LOW PROTEIN rats (burst number) - casein vs. maltodextrin')
print(ro.r('lp_casvmalt'))

# Analysis of protein preference

ro.globalenv['nppref'] = df2[df2['diet'] == 'np']
ro.globalenv['lppref'] = df2[df2['diet'] != 'np']

ro.r('proteinPref = t.test(nppref[\'pref\'], lppref[\'pref\'], paired=FALSE)')
print(ro.r('proteinPref'))


#from rpy2.robjects.packages import importr
#from rpy2.robjects.vectors import StrVector
#
#stats = importr('stats')
#base = importr('base')
#datasets = importr('datasets')

#myresult = stats.t_test(data_pref1['cas_pref'][msk], data_pref1['cas_pref'][~msk],
#                        **{'paired': False})
#

#r_df = pandas2ri.py2ri(df[['ratid', 'sol', 'diet', 'total', 'bMean', 'bNum']])
#ro.r('r_df[, \'sol\'] <- as.factor(r_df[, \'sol\'])')
#ro.r('r_df[, \'diet\'] <- as.factor(r_df[, \'diet\'])')
#ro.r('r_df[, \'ratid\'] <- as.factor(r_df[, \'ratid\'])')

#ro.globalenv['r_df'] = ro.r(pandas2ri.py2ri(df[['ratid', 'sol', 'diet', 'total', 'bMean', 'bNum']]))
#
#
#ro.r('result = aov(formula = licks ~ sol * diet + Error(ratid / sol), data = cas9_pref1_r)')
#myresult2 = ro.r('summary(result)')
#
#print(myresult2)
#
##myresult2 = stats.aov(data ~ sol * diet)
#
##z = pd.DataFrame([vars(x) for x in rats.preference1_cas])

#filename = userhome + '\\Dropbox\\Python\\cas9\\cas9_stats\\cas9_pref1_r2.csv'
#r_df.to_csv(filename)
#ro.r('library(readr)')
#ro.r('cas9_pref1_r <- read_csv("C:/Users/James Rig/Dropbox/Python/cas9/cas9_stats/cas9_pref1_r2.csv",col_types = cols(diet = col_factor(levels = c("np","lp")), sol = col_factor(levels = c("c","m"))))')
#