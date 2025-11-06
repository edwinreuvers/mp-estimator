# -*- coding: utf-8 -*-
"""
 Example 4:
     Estimate force-length and CE force-velocity parameters using the 
     'improved method'

 @author: Edwin D.H.M. Reuvers
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
    optsSR = {
        'dataDir':          dataDirMus+'SR',
        'iCols':            (0,1,3),
        'idxSRcon':         'auto',
        'nSRsamp':          'auto',
        'dispFig':          False,
        }
    
    #%% Load 'default' parameter values and data
    defpar = loaddata.get_muspar()
    dataQR,idxQRmin,idxQRpre,idxQRpst = loaddata.qr(optsQR)
    dataSR,idxSRcon = loaddata.sr(optsSR)
    
    #%% Estimate contraction dynamics parameter values using the 'improved method'
    estpar,dataQRout,dataSRout = estimate.im(dataQR,dataSR,defpar,do_print=True)
    filepath = os.path.join(parDir,mus,f'{mus}_example4.pkl')
    pickle.dump(estpar, open(filepath, 'wb'))