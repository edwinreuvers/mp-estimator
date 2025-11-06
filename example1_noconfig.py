# -*- coding: utf-8 -*-
"""
 Example 1:
     Estimate CE, SEE and PEE force-length parameters only.
     
     This version has no pre-defined configurations settings. To use the script:
         1. Select datafiles: 
             In the first pop-up window, choose the CSV file containing your data.
        2. Select data columns:
            In the second pop-up window, specify which columns represent:
                Time (time points for the data)
                Muscle Length (muscle length at each time point)
                Muscle Force (force exerted by the muscle at each time point)
        3. Review data selection:
            In the third pop-up window, check the data selection to ensure the columns and file are correct.
"""

#%% Load packages
from mp_estimator import estimate, loaddata
 
#%% Load 'default' parameter values and data
defpar = loaddata.get_muspar()
dataQR,idxQRmin,idxQRpre,idxQRpst = loaddata.qr()

#%% Estimate parameter values of CE, SEE and PEE force-length parameters
estpar,dataQRout = estimate.fl(dataQR,defpar)