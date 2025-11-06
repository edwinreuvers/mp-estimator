# -*- coding: utf-8 -*-
"""
 Example 2:
     Estimate CE force-velocity parameters only
"""

#%% Load packages
import os, pickle

# Set directories
cwd = os.path.dirname(os.path.abspath(__file__))
baseDir = os.path.join(cwd,'')
dataDir = os.path.join(baseDir,'data')
parDir = os.path.join(baseDir,'parameters')

# Load custom functions
from mp_estimator import estimate, loaddata

#%% Muscle & folders
for mus in ['GMs1', 'GMs2', 'GMs3']:
    dataDirMus = os.path.join(dataDir,mus,'')
    
    #%% Settings
    optsSR = {
        'dataDir':          dataDirMus+'SR',
        'iCols':            (0,1,3),
        'idxSRcon':         'auto',
        'nSRsamp':          'auto',
        'dispFig':          False,
        }
    
    #%% Load 'default' parameter values and data
    defpar = loaddata.get_muspar()
    dataSR,idxSRcon = loaddata.sr(optsSR)
    
    #%% Estimate muscle parameter values of CE force-velocity parameters
    estpar,dataSRout = estimate.fv(dataSR,defpar)
    filepath = os.path.join(parDir,mus,f'{mus}_example2.pkl')
    pickle.dump(estpar, open(filepath, 'wb'))