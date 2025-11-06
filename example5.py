# -*- coding: utf-8 -*-
"""
 Example 5:
     Estimate (de-)activation dynamics parameters

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
# for mus in ['GMs1', 'GMs2', 'GMs3']:
for mus in ['GMe1']:
    dataDirMus = os.path.join(dataDir,mus,'')
    
    #%% Settings
    optsISOM = {
        'dataDir':          dataDirMus+'ISOM',
        'iCols':            (0,1,3,2),
        'durStimOffset':    0.1,
        'dispFig':          False,
        }
    
    #%% Load 'default' parameter values and data
    defpar = loaddata.get_muspar()
    dataISOM, idxSEL = loaddata.isom(optsISOM)
    
    #%% Load parameter values from 'improved method'
    filepath = os.path.join(parDir,mus,f'{mus}_example4.pkl')
    estpar = pickle.load(open(filepath, 'rb'))
    
    #%% Estimate excitation dynamics parameter values
    estpar,dataACTout = estimate.act(dataISOM,defpar,estpar,do_print=True)
    filepath = os.path.join(parDir,mus,f'{mus}_example5_IM.pkl')
    pickle.dump(estpar, open(filepath, 'wb'))