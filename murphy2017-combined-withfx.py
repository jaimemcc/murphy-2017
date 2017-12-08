# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:32:27 2017

@author: LocalAdmin1
"""

# Uncomment these imports for R statistics
makefigs = True
savefigs = True
statson = True

if statson == True:
    import rpy2.robjects as ro
    from rpy2.robjects import r, pandas2ri, numpy2ri
    pandas2ri.activate()
    numpy2ri.activate()
    from scipy import stats

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

import JM_custom_figs as jmfig

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
        licks = medfilereader(self.medfile,
                                  varsToExtract = sub2var(self, substance),
                                                    remove_var_header = True)
        lickData = lickCalc(licks, burstThreshold=0.5, binsize=120)        
        
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
                            
def metafilereader(filename):
    
    f = open(filename, 'r')
    f.seek(0)
    header = f.readlines()[0]
    f.seek(0)
    filerows = f.readlines()[1:]
    
    tablerows = []
    
    for i in filerows:
        tablerows.append(i.split('\t'))
        
    header = header.split('\t')
    # need to find a way to strip end of line \n from last column - work-around is to add extra dummy column at end of metafile
    return tablerows, header

def isnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')

def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [isnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn

"""
This function will calculate data for bursts from a train of licks. The threshold
for bursts and clusters can be set. It returns all data as a dictionary.
"""
def lickCalc(licks, offset = [], burstThreshold = 0.25, runThreshold = 10, 
             binsize=60, histDensity = False):
    # makes dictionary of data relating to licks and bursts
    if type(licks) != np.ndarray or type(offset) != np.ndarray:
        try:
            licks = np.array(licks)
            offset = np.array(offset)
        except:
            print('Licks and offsets need to be arrays and unable to easily convert.')
            return

    lickData = {}
    
    if len(offset) > 0:
        lickData['licklength'] = offset - licks
        lickData['longlicks'] = [x for x in lickData['licklength'] if x > 0.3]
    else:
        lickData['licklength'] = []
        lickData['longlicks'] = []

    lickData['licks'] = np.concatenate([[0], licks])
    lickData['ilis'] = np.diff(lickData['licks'])
    lickData['freq'] = 1/np.mean([x for x in lickData['ilis'] if x < burstThreshold])
    lickData['total'] = len(licks)
    
    # Calculates start, end, number of licks and time for each BURST 
    lickData['bStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]  
    lickData['bInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)]
    lickData['bEnd'] = [lickData['licks'][i-1] for i in lickData['bInd'][1:]]
    lickData['bEnd'].append(lickData['licks'][-1])

    lickData['bLicks'] = np.diff(lickData['bInd'] + [len(lickData['licks'])])    
    lickData['bTime'] = np.subtract(lickData['bEnd'], lickData['bStart'])
    lickData['bNum'] = len(lickData['bStart'])
    if lickData['bNum'] > 0:
        lickData['bMean'] = np.nanmean(lickData['bLicks'])
    else:
        lickData['bMean'] = 0
    
    lickData['bILIs'] = [x for x in lickData['ilis'] if x > burstThreshold]

    # Calculates start, end, number of licks and time for each RUN
    lickData['rStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]  
    lickData['rInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > runThreshold)]
    lickData['rEnd'] = [lickData['licks'][i-1] for i in lickData['rInd'][1:]]
    lickData['rEnd'].append(lickData['licks'][-1])

    lickData['rLicks'] = np.diff(lickData['rInd'] + [len(lickData['licks'])])    
    lickData['rTime'] = np.subtract(lickData['rEnd'], lickData['rStart'])
    lickData['rNum'] = len(lickData['rStart'])

    lickData['rILIs'] = [x for x in lickData['ilis'] if x > runThreshold]
    try:
        lickData['hist'] = np.histogram(lickData['licks'][1:], 
                                    range=(0, 3600), bins=int((3600/binsize)),
                                          density=histDensity)[0]
    except TypeError:
        print('Problem making histograms of lick data')
        
    return lickData

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

def condhistFig(ax, df, factor, sol='maltodextrin'):
    if sol == 'casein':
        NRcolor = 'black'
        PRcolor = 'xkcd:kelly green'
    else:
        NRcolor = 'xkcd:bluish grey'
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
                 grouplabel=['Day 1', 'Day 2'],
                 scattersize = 30,
                 ax=ax)


metafile = userhome + '\\Documents\\GitHub\\murphy-2017\\CAS9_metafile.txt'
metafileData, metafileHeader = metafilereader(metafile)

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

## Body weight and food intake

bwmetafile = userhome + '\\Documents\\GitHub\\murphy-2017\\CAS9bw_metafile.csv'

data = pd.read_csv(bwmetafile, index_col=['rat', 'diet'])

# Statistics

if statson == True:
    data = data[:].stack()
    data = data.to_frame()
    data.reset_index(inplace=True) 
    data.columns = ['rat', 'diet', 'day', 'licks']
    ro.globalenv['r_data'] = data
    
    ro.r('bodyweight = aov(formula = licks ~ day * diet + Error(rat / day), data = r_data)')
    
    print(ro.r('summary(bodyweight)'))


data = pd.read_csv(bwmetafile, index_col=['rat'])

np_mean = data[data['diet'] == 'np'].mean()
np_sem = data[data['diet'] == 'np'].std() / np.sqrt(len(data['diet'] == 'np'))

lp_mean = data[data['diet'] == 'lp'].mean()
lp_sem = data[data['diet'] == 'lp'].std() / np.sqrt(len(data['diet'] == 'lp'))

# Figure 1A - Body weight
if makefigs == True:
    fig1a = plt.figure(figsize=(3.2,2.4))
    ax = plt.subplot(1,1,1)
    np_mean.plot(yerr=np_sem, color='xkcd:charcoal', marker='o', markerfacecolor='white')
    lp_mean.plot(yerr=lp_sem, color='xkcd:kelly green', marker='o', markerfacecolor='white')
    ax.set_ylim([400, 550])
    ax.set_xlim([-1, 17])
    plt.xticks([1,6,11,16], ('0', '5', '10', '15'))
    plt.yticks([400, 450, 500, 550])
    ax.set_ylabel('Body weight (g)')
    ax.set_xlabel('Days since diet switch')

# Figure 1B - Food intake
# Data from CAS9_exptdetails.xls > fi_cas9+cas56

foodintake_np = [22.0, 25.6, 27.6, 26.1]
foodintake_lp = [25.5, 27.3, 29.1, 28.5]

foodintake_npAll = [22.0, 25.6, 27.6, 26.1, 22.9, 21.0, 21.3, 20.6]
foodintake_lpAll = [25.5, 27.3, 29.1, 28.5, 27.4, 29.9, 23.5, 26.3]

if makefigs == True:
    
    fi = data2obj1D([foodintake_np, foodintake_lp])
    mpl.rcParams['figure.subplot.left'] = 0.25
    fig1b = plt.figure(figsize=(1.8,2.4))
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

    fi = data2obj1D([foodintake_npAll, foodintake_lpAll])
    fig1b2 = plt.figure(figsize=(1.8,2.4))
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


mpl.rcParams['figure.subplot.left'] = 0.15

if statson == True:
    fi_stats = stats.ttest_ind(foodintake_npAll, foodintake_lpAll)
    print(fi_stats)

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

df3 = df[['ratid', 'diet', 'day', 'total']]

### Need to append preference data to this df so that I can compare different
### amounts of licks

if makefigs == True:
    
    mpl.rcParams['figure.subplot.wspace'] = 0.1
    mpl.rcParams['figure.subplot.left'] = 0.15
    fig2a, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(3.2, 2.4))
    
    condhistFig(ax[0], dfc1, 'hist', sol='casein')
    fig2a.text(0.55, 0.04, 'Time (minutes)', ha='center')
    ax[0].set_ylabel('Licks per 2 min')
    ax[0].set_ylim([-20, 410])
    
    condhistFig(ax[1], dfc2, 'hist', sol='casein')
    
    fig2b, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(3.2, 2.4))
    cond2Dfig(ax, df, 'total', sol='casein')    
    plt.yticks([0, 2000, 4000, 6000], ('0', '2', '4', '6'))
    ax.set_ylim([0, 7800])
    ax.set_ylabel('Licks (x1000)')
    
    fig2c, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(3.2, 2.4))
    condhistFig(ax[0], dfm1, 'hist')
    fig2c.text(0.55, 0.04, 'Time (minutes)', ha='center')
    ax[0].set_ylabel('Licks per 2 min')
    ax[0].set_ylim([-20, 410])
    
    condhistFig(ax[1], dfm2, 'hist')

    fig2d, ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True, figsize=(3.2, 2.4))
    cond2Dfig(ax, df, 'total')
    plt.yticks([0, 2000, 4000, 6000], ('0', '2', '4', '6'))
    ax.set_ylim([0, 7800])
    ax.set_ylabel('Licks (x1000)')

# Figure of total licks during conditioning
    fig2e = plt.figure(figsize=(3.2, 2.4))
    ax = plt.subplot(1,1,1)
    nplp2Dfig(df2, 'total', ax=ax)
    plt.yticks([0, 5000, 10000, 15000], ('0', '5', '10', '15'))
    ax.set_ylabel('Licks (x1000)')

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
    
    ro.r('nrpr_cas_DAY1 = t.test(r_df$total[r_df$diet=="lp" & r_df$cday=="1"], r_df$total[r_df$diet=="np" & r_df$cday=="1"], paired=FALSE)')
    print('LOW vs NORMAL PROTEIN rats (licks per CASEIN conditioning) - DAY 1')
    print(ro.r('nrpr_cas_DAY1'))
    
    ro.r('nrpr_cas_DAY2 = t.test(r_df$total[r_df$diet=="lp" & r_df$cday=="2"], r_df$total[r_df$diet=="np" & r_df$cday=="2"], paired=FALSE)')
    print('LOW vs NORMAL PROTEIN rats (licks per CASEIN conditioning) - DAY 2')
    print(ro.r('nrpr_cas_DAY2'))
    
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

# Combining dfs from conditioning and preference day to look at changes in total licks
df4 = df[['ratid', 'diet']][:24]
df4.insert(2, 'cday', [5]*24)
prefdaytotal = (df.total[:24]+df.total[24:])
df4.insert(3,'total',prefdaytotal)

df5 = pd.concat([df3, df4])

#if statson == True:
#daymsk = df.day == 'c'
#r_df = df5[['ratid', 'diet', 'day', 'total']][solmsk]
#ro.globalenv['r_df'] = r_df
#    
#    ro.r('casein_cond = aov(formula = total ~ cday * diet + Error(ratid / cday), data = r_df)')
#    print('Casein during conditioning')
#    print(ro.r('summary(casein_cond)'))

# Figure 3A - Licks over time, histogram
if makefigs == True:
    mpl.rcParams['figure.subplot.wspace'] = 0.1
    mpl.rcParams['figure.subplot.left'] = 0.15
    fig3a, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(3.2, 2.4))
    
    prefhistFig(ax[0], ax[1], df, 'hist')
    fig3a.text(0.55, 0.04, 'Time (min)', ha='center')
    ax[0].set_ylabel('Licks per 2 min')

# Figure 3B - casein licks vs maltodextrin licks
    mpl.rcParams['figure.subplot.left'] = 0.25
    fig3b = plt.figure(figsize=(2.4, 2.4))
    ax = plt.subplot(1,1,1)                
    casVmaltFig(ax, df)
    ax.set_xlabel('Licks for casein')
    ax.set_ylabel('Licks for maltodextrin')
    plt.yticks([0, 2000, 4000, 6000])
       
    mpl.rcParams['figure.subplot.left'] = 0.15

# Figure 3C - casein and malt preference (licks)
# Analysis of licks
    fig3c = plt.figure(figsize=(3.2, 2.4))
    ax = plt.subplot(1,1,1) 
    nplp2Dfig(df, 'total', ax)
    ax.set_ylabel('Licks (x1000)')
    plt.yticks([0, 2000, 4000, 6000], ('0', '2', '4', '6'))

# Analysis of palatability

# Figure 4A - licks per burst

    fig4a = plt.figure(figsize=(3.2, 2.4))
    ax = plt.subplot(1,1,1)
    nplp2Dfig(df, 'bMean', ax)
    ax.set_ylabel('Average licks per cluster')
    ax.set_yticks([0, 50, 100, 150])
    
    fig4b = plt.figure(figsize=(3.2, 2.4))
    ax = plt.subplot(1,1,1)
    nplp2Dfig(df, 'bNum', ax)
    ax.set_ylabel('Number of clusters')
    ax.set_yticks([0, 50, 100, 150, 200])

# Figure 3D - protein preference

dietmsk = df2.diet == 'np'
a = data2obj1D([df2['pref'][dietmsk], df2['pref'][~dietmsk]])

if makefigs == True:
    mpl.rcParams['figure.subplot.left'] = 0.25
    fig3d = plt.figure(figsize=(1.8, 2.4))
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
    ro.r('proteinPref')


if savefigs == True:
    savefolder = userhome + '\\Desktop\\Temp Figs\\'
    
    fig1a.savefig(savefolder + '01a_bodyweight.eps')
    fig1b.savefig(savefolder + '01b_foodintake.eps')
    fig1b2.savefig(savefolder + '01b2_foodintake.eps')
    
    fig2a.savefig(savefolder + '02a_condhist_cas.eps')
    fig2b.savefig(savefolder + '02b_condbar_cas.eps')
    fig2c.savefig(savefolder + '02c_condhist_malt.eps')
    fig2d.savefig(savefolder + '02d_condbar_malt.eps')
    fig2e.savefig(savefolder + '02e_condtotallicks.eps')
       
    fig3a.savefig(savefolder + '03a_prefhist.eps')
    fig3b.savefig(savefolder + '03b_casvmalt.eps')
    fig3c.savefig(savefolder + '03c_preftotal.eps')
    fig3d.savefig(savefolder + '03d_caseinpref.eps')
    
    fig4a.savefig(savefolder + '04a_prefburstmean.eps')
    fig4b.savefig(savefolder + '04b_prefburstnum.eps')

