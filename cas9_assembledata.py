# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 10:32:27 2017

@author: LocalAdmin1
"""
import sys
sys.path.insert(0,'C:\\Github\\functions-and-figures\\')

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

import JM_general_functions as jmf
import JM_custom_figs as jmfig

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

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
#import rpy2.robjects as ro

#from rpy2.robjects import r, pandas2ri, numpy2ri
#pandas2ri.activate()
#numpy2ri.activate()

import os
import timeit

tic = timeit.default_timer()

userhome = os.path.expanduser('~')
datafolder = 'data\\'
        
class Session(object):
    
    def __init__(self, sessionID, metafiledata, hrows, datafolder):
        self.hrow = hrows
        self.sessionID = sessionID
        self.rat = metafiledata[hrows['rat']].replace('.', '-')
        self.session = metafiledata[hrows['session']]
        self.medfile = datafolder + metafiledata[hrows['medfile']]
        self.sessioncode = self.hrow['sessioncode']

        self.bottleA = self.hrow['bottleA']
        self.bottleB = self.hrow['bottleB']
                    
    def extractlicks(self, substance):
        licks = jmf.medfilereader(self.medfile,
                                  varsToExtract = sub2var(self, substance),
                                                    remove_var_header = True)
        lickData = jmf.lickCalc(licks, burstThreshold=0.5, minburstlength=3, binsize=120)        
        
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

def metafile2sessions(metafile, datafolder):
#    jmf.metafilemaker(xlfile, metafile, sheetname=sheetname, fileformat='txt')
    rows, header = jmf.metafilereader(metafile)
    
    hrows = {}
    for idx, field in enumerate(header):
        hrows[field] = idx
    
    sessions = {}
    
    for row in rows:
        sessionID = row[hrows['rat']].replace('.','-') + '_' + row[hrows['session']]
        sessions[sessionID] = Session(sessionID, row, hrows, datafolder)
    
    return sessions   

for session in sessions:
      
    x = sessions[session]
    x.extractlicks('casein')

metafile = 'CAS9_metafile.txt'
#metafileData, metafileHeader = jmf.metafilereader(metafile)

sessions = metafile2sessions(metafile, datafolder)

#exptsuffix = ''
#includecol = 10
#
#try:
#    type(rats2)
#    print('Using existing data')
#except NameError:
#    print('Assembling data from Med Associates files')
#    
#    rats = {}
#    
#    for i in metafileData:
#        if int(i[includecol]) == 1:
#            rowrat = str(i[1])
#            if rowrat not in rats:
#                rats[rowrat] = Rat(rowrat)
#            rats[rowrat].loadsession(i, metafileHeader)
#            
#    for i in rats:
#        for j in rats[i].sessions:
#    #        print('Analysing rat ' + i + ' in session ' + j)
#            x = rats[i].sessions[j]
#            try:
#                x.lickData_cas = x.extractlicks('casein')
#            except IndexError:
#                print('Difficulty extracting casein licks')
#                
#            try:
#                x.lickData_malt = x.extractlicks('maltodextrin')
#            except IndexError:
#                print('Difficulty extracting casein licks')
#                
#                
#            x.lickData_sacc = x.extractlicks('saccharin')
#            
#            x.designatesession()
