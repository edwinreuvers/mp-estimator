# -*- coding: utf-8 -*-
"""
 Example 1:
     Estimate CE, SEE and PEE force-length parameters only.
     
     This version has pre-defined configuration settings. 
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
    optsQR = {
        'dataDir':          dataDirMus+'QR',
        'iCols':            (0,1,3),
        'idxQRmin':         'auto',
        'idxQRpre':         'auto',
        'idxQRpst':         'auto',
        'nQRsamp':          'auto',
        'dispFig':          False,
        }
        
    #%% Load 'default' parameter values and data
    defpar = loaddata.get_muspar()
    dataQR,idxQRmin,idxQRpre,idxQRpst = loaddata.qr(optsQR)
    
    #%% Estimate parameter values of CE, SEE and PEE force-length parameters
    defpar = loaddata.get_muspar()
    estpar,dataQRout = estimate.fl(dataQR,defpar)
    filepath = os.path.join(parDir,mus,f'{mus}_example1.pkl')
    pickle.dump(estpar, open(filepath, 'wb'))