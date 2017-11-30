# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:32:27 2017

@author: LocalAdmin1
"""

# Uncomment these imports for R statistics
makefigs = True
savefigs = False
statson = False

if statson == True:
    import rpy2.robjects as ro
    from rpy2.robjects import r, pandas2ri, numpy2ri
    pandas2ri.activate()
    numpy2ri.activate()

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

import JM_general_functions as jmf
import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

#plt.style.use('murphy')

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (3.2, 2.4)
mpl.rcParams['figure.dpi'] = 100

mpl.rcParams['font.size'] = 8.0
mpl.rcParams['axes.labelsize'] = 'medium'
mpl.rcParams['ytick.labelsize'] = 'small'

mpl.rcParams['figure.subplot.left'] = 0.15
mpl.rcParams['figure.subplot.bottom'] = 0.20

mpl.rcParams['errorbar.capsize'] = 5

mpl.rcParams['savefig.transparent'] = True

mpl.rcParams['axes.spines.top']=False
mpl.rcParams['axes.spines.right']=False

## Colour scheme
col={}
col['np_cas'] = 'xkcd:silver'
col['np_malt'] = 'white'
col['lp_cas'] = 'xkcd:kelly green'
col['lp_malt'] = 'xkcd:light green'

import pandas as pd

import os
import timeit

tic = timeit.default_timer()

userhome = os.path.expanduser('~')
datafolder = userhome + '\\Documents\\GitHub\\murphy-2017\\cas9_medfiles\\'

class Rat(object):
    
    nRats = 0
    nSessions = 0
    
    def __init__(self, data):      
        self.rat = data
        self.sessions = {}
        
        Rat.nRats += 1
                
    def loadsession(self, data, header):
        self.session = 's'+str(data[2]) #should reference column of data with session number
        self.sessions[self.session] = Session(data, header, self.rat, self.session)
       
        Rat.nSessions += 1
        
class Session(object):
    
    def __init__(self, data, header, rat, session):
        self.hrow = {}
        for idx, col in enumerate(header):
            self.hrow[col] = data[idx]
        self.medfile = datafolder + self.hrow['medfile']
        self.sessioncode = self.hrow['sessioncode']
        self.rat = str(rat)
        self.session = session
        self.bottleA = self.hrow['bottleA']
        self.bottleB = self.hrow['bottleB']
        
        if hasattr(rats[self.rat], 'diet'):
            if rats[self.rat].diet != self.hrow['diet']:
                print('Wrong diet for rat, must be a mistake in metafile')
        else:
            rats[self.rat].diet = self.hrow['diet']
                    
    def extractlicks(self, substance):
        licks = jmf.medfilereader(self.medfile,
                                  varsToExtract = sub2var(self, substance),
                                                    remove_var_header = True)
        lickData = jmf.lickCalc(licks, burstThreshold=0.5, binsize=120)        
        
        return lickData

    def designatesession(self):
        if self.sessioncode == 'casein1':
            if hasattr(rats[self.rat], 'casein1'):
                print('Casein 1 data already added. Check metafile for duplication.')
            else:
                rats[self.rat].casein1 = self.lickData_cas
                    
        if self.sessioncode == 'casein2':
            if hasattr(rats[self.rat], 'casein2'):
                print('Casein 2 data already added. Check metafile for duplication.')
            else:
                rats[self.rat].casein2 = self.lickData_cas

        if self.sessioncode == 'maltodextrin1':
            if hasattr(rats[self.rat], 'maltodextrin1'):
                print('Maltodextrin 1 data already added. Check metafile for duplication.')
            else:
                rats[self.rat].maltodextrin1 = self.lickData_malt
                    
        if self.sessioncode == 'maltodextrin2':
            if hasattr(rats[self.rat], 'maltodextrin2'):
                print('Maltodextrin 2 data already added. Check metafile for duplication.')
            else:
                rats[self.rat].maltodextrin2 = self.lickData_malt
                    
        if self.sessioncode == 'preference1':
            if hasattr(rats[self.rat], 'preference1-cas'):
                print('Preference 1 data already added. Check metafile for duplication.')
            else:
                rats[self.rat].preference1_cas = self.lickData_cas
                rats[self.rat].preference1_malt = self.lickData_malt
                            
def sub2var(session, substance):
    varsOut = []
    if session.bottleA == substance:
        varsOut.append('b')
    if session.bottleB == substance:
        varsOut.append('d')
    return varsOut

def data2obj1D(data):
    obj = np.empty(len(data), dtype=np.object)
    for i,x in enumerate(data):
        obj[i] = np.array(x)  
    return obj

def data2obj2D(data):
    obj = np.empty((np.shape(data)[0], np.shape(data)[1]), dtype=np.object)
    for i,x in enumerate(data):
        for j,y in enumerate(x):
            obj[i][j] = np.array(y)
    return obj

def nplp2Dfig(df, key, ax):
    dietmsk = df.diet == 'np'
    casmsk = df.sol == 'c'
    
    a = [[df[key][dietmsk & casmsk], df[key][dietmsk & ~casmsk]],
          [df[key][~dietmsk & casmsk], df[key][~dietmsk & ~casmsk]]]

    x = data2obj2D(a)

    ax, x, _, _ = jmfig.barscatter(x, paired=True,
                 barfacecoloroption = 'individual',
                 barfacecolor = [col['np_cas'], col['np_malt'], col['lp_cas'], col['lp_malt']],
                 scatteredgecolor = ['xkcd:charcoal'],
                 scatterlinecolor = 'xkcd:charcoal',
                 grouplabel=['NR', 'PR'],
                 scattersize = 60,
                 ax=ax)

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

metafile = userhome + '\\Documents\\GitHub\\murphy-2017\\CAS9_metafile.txt'
metafileData, metafileHeader = jmf.metafilereader(metafile)

exptsuffix = ''
includecol = 10

try:
    type(rats)
    print('Using existing data')
except NameError:
    print('Assembling data from Med Associates files')
    
    rats = {}
    
    for i in metafileData:
        if int(i[includecol]) == 1:
            rowrat = str(i[1])
            if rowrat not in rats:
                rats[rowrat] = Rat(rowrat)
            rats[rowrat].loadsession(i, metafileHeader)
            
    for i in rats:
        for j in rats[i].sessions:
    #        print('Analysing rat ' + i + ' in session ' + j)
            x = rats[i].sessions[j]
            
            x.lickData_cas = x.extractlicks('casein')
            x.lickData_malt = x.extractlicks('maltodextrin')
            x.lickData_sacc = x.extractlicks('saccharin')
            
            x.designatesession()

## Analysis of conditioning days
# Assembling conditioning data

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
df.insert(4,'cday',[1]*24 + [2]*24 + [1]*24 + [2]*24)

df2 = df[['ratid', 'diet']][:48]
casall = df[df['day'] == 'c1']['total'] + df[df['day'] == 'c2']['total']
maltall = df[df['day'] == 'm1']['total'] + df[df['day'] == 'm2']['total']

df2.insert(2,'sol',['c']*24 + ['m']*24)
df2.insert(3,'total', casall.append(maltall))

df3 = df[['ratid', 'diet', 'sol', 'total']]

### Need to append preference data to this df so that I can compare different
### amounts of licks

def condhistFig(ax, df, factor, sol='maltodextrin'):
    if sol == 'casein':
        NRcolor = 'black'
        PRcolor = 'xkcd:kelly green'
    else:
        NRcolor = 'xkcd:silver'
        PRcolor = 'xkcd:light green'
        
    dietmsk = df.diet == 'np'

    shadedError(ax, df[factor][dietmsk], linecolor=NRcolor)
    ax = shadedError(ax, df[factor][~dietmsk], linecolor=PRcolor)
    ax.set_xticks([0,10,20,30])
    ax.set_xticklabels(['0', '20', '40', '60'])
    
def cond2Dfig(ax, df, factor, sol='maltodextrin'):
    if sol == 'casein':
        day1msk = df.day == 'c1'
        day2msk = df.day == 'c2'
    else:
        day1msk = df.day == 'm1'
        day2msk = df.day == 'm2'
        
    dietmsk = df.diet == 'np'
   
    a = [[df[factor][day1msk & dietmsk], df[factor][day1msk & ~dietmsk]],
          [df[factor][day2msk & dietmsk], df[factor][day2msk & ~dietmsk]]]

    x = data2obj2D(a)
    
    if sol == 'casein':
        barfacecolor = [col['np_cas'], col['lp_cas'], col['np_cas'], col['lp_cas']]
    else:
        barfacecolor = [col['np_malt'], col['lp_malt'], col['np_malt'], col['lp_malt']]
        
    ax, x, _, _ = jmfig.barscatter(x, paired=False,
                 barfacecoloroption = 'individual',
                 barfacecolor = barfacecolor,
                 scatteredgecolor = ['xkcd:charcoal'],
                 scatterlinecolor = 'xkcd:charcoal',
                 grouplabel=['NR', 'PR'],
                 scattersize = 60,
                 ax=ax)
    
mpl.rcParams['figure.subplot.wspace'] = 0.1
mpl.rcParams['figure.subplot.left'] = 0.15
fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(4.8, 3.2))

condhistFig(ax[0], dfc1, 'hist', sol='casein')
fig.text(0.55, 0.04, 'Time (minutes)', ha='center')
ax[0].set_ylabel('Licks per 2 min')

condhistFig(ax[1], dfc2, 'hist', sol='casein')
fig.text(0.55, 0.04, 'Time (minutes)', ha='center')

fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(3.2, 3.2))
cond2Dfig(ax, df, 'total', sol='casein')

fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(4.8, 3.2))
condhistFig(ax[0], dfm1, 'hist')
fig.text(0.55, 0.04, 'Time (minutes)', ha='center')
ax[0].set_ylabel('Licks per 2 min')

condhistFig(ax[1], dfm2, 'hist')
fig.text(0.55, 0.04, 'Time (minutes)', ha='center')

fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(3.2, 3.2))
cond2Dfig(ax, df, 'total')

## Statistics with R - conditioning
# Day 1 vs 2, PR vs NR for CASEIN

if statson == True:
    solmsk = df.sol == 'c'
    r_df = df[['ratid', 'diet', 'cday', 'total']][solmsk]
    ro.globalenv['r_df'] = r_df
    
    ro.r('casein_cond = aov(formula = total ~ cday * diet + Error(ratid / cday), data = r_df)')
    print('Casein during conditioning')
    print(ro.r('summary(casein_cond)'))
    
    ro.r('pr_cas_day12 = t.test(r_df$total[r_df$diet=="lp" & r_df$cday=="1"], r_df$total[r_df$diet=="lp" & r_df$cday=="2"], paired=TRUE)')
    print('LOW PROTEIN rats (licks per CASEIN conditioning) - day 1 vs. day 2')
    print(ro.r('pr_cas_day12'))
    
    ro.r('nr_cas_day12 = t.test(r_df$total[r_df$diet=="np" & r_df$cday=="1"], r_df$total[r_df$diet=="np" & r_df$cday=="2"], paired=TRUE)')
    print('NORMAL PROTEIN rats (licks per CASEIN conditioning) - day 1 vs. day 2')
    print(ro.r('nr_cas_day12'))
    
    # Day 1 vs 2, PR vs NR for CASEIN
    solmsk = df.sol == 'm'
    r_df = df[['ratid', 'diet', 'cday', 'total']][solmsk]
    ro.globalenv['r_df'] = r_df
    
    ro.r('malto_cond = aov(formula = total ~ cday * diet + Error(ratid / cday), data = r_df)')
    print('Maltodextrin during conditioning')
    print(ro.r('summary(malto_cond)'))
    
    #ro.r('pr_malt_day12 = t.test(r_df$total[r_df$diet=="lp" & r_df$cday=="1"], r_df$total[r_df$diet=="lp" & r_df$cday=="2"], paired=TRUE)')
    #print('LOW PROTEIN rats (licks per MALTO conditioning) - day 1 vs. day 2')
    #print(ro.r('pr_malt_day12'))
    #
    #ro.r('nr_malt_day12 = t.test(r_df$total[r_df$diet=="np" & r_df$cday=="1"], r_df$total[r_df$diet=="np" & r_df$cday=="2"], paired=TRUE)')
    #print('NORMAL PROTEIN rats (licks per MALTO conditioning) - day 1 vs. day 2')
    #print(ro.r('nr_malt_day12'))

# Figure of total licks during conditioning

if makefigs == True:
    fig = plt.figure(figsize=(3.2, 2.4))
    ax = plt.subplot(1,1,1)
    nplp2Dfig(df2, 'total', ax=ax)
    plt.yticks([0, 5000, 10000, 15000], ('0', '5', '10', '15'))
    ax.set_ylabel('Licks (x1000)')
    #plt.savefig('ADD FILEPATH/03_condlicks.eps')
    plt.title('Consumption during conditioning')

if statson == True:
    r_df = df2[['ratid', 'sol', 'diet', 'total']]
    ro.globalenv['r_df'] = r_df
    ro.r('condlicks = aov(formula = total ~ sol * diet + Error(ratid / sol), data = r_df)')
    print('Licks during conditioning')
    print(ro.r('summary(condlicks)'))

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
if makefigs == True:
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

if makefigs == True:
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
if statson == True:
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

