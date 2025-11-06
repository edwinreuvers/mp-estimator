# -*- coding: utf-8 -*-
"""
    loaddata.py contains all functions to load data of the QR, SR and ISOM
    experiments.

    author: Edwin D.H.M. Reuvers
    v1.0.0, september 2025
"""

#%% Import packages
import glob, os
import numpy as np
import pandas as pd
from scipy import signal
from .interface import DataPlot, OpenFileDialog, SelDataCol
from .helpers import get_stim

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_muspar():
    import numpy as np
    muspar = {}

    # Force-length relationship
    muspar['w'] = 0.5
    muspar['n'] = 2

    # Force-velocity relationship
    muspar['fasymp'] = 1.5
    muspar['slopfac'] = 2
    muspar['vfactmin'] = 0.1

    # (De-)activation dynamics
    muspar['q0'] = 0.005
    muspar['gamma_0'] = 1e-5
    muspar['kCa'] = 8e-6
    muspar['tact'] = 20e-3
    muspar['tdeact'] = 20e-3
    muspar['a_act'] = -7.369163055099003
    muspar['b_act'] = np.array([5.170927028993413, 0.5955111970420514, 0])
    
    return muspar

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def qr(opts={},fQR={}):
    if fQR == {}:
        fQR,opts = get_qr_opts(opts)
    
    if len(fQR)<3:
        raise ValueError("Less than three datafiles are selected, at least 3 datafiles are necessary.") 
        return None, None, None, None
    
    iCols = opts['iCols']
    nFiles = opts['nFiles']
    nQRsamp = opts['nQRsamp']
    
    # Make listed-dic to store data   
    dataQR = [{'time': None, 'lmtc': None, 'fsee': None, 'filename': None, 'idxQRmin': None, 'idxQRpre': None, 'idxQRpst': None} for _ in range(nFiles)]
    # Now select data etc
    for iFile,filepath in enumerate(fQR):
        filename = filepath.rsplit('/', 1)[-1][:-4]
        data = pd.read_csv(filepath).T.to_numpy()
        time = data[iCols[0]]
        lmtc = data[iCols[1]] 
        fsee = data[iCols[2]]
        
        # First moving average filter (50 ms)
        dt = np.diff(time).mean()
        N = int(np.round(1/(20*dt)))
        B = np.ones(N)/N
        fseeFilt = signal.filtfilt(B,1,fsee)

        # Find indices where SEE force is minimum
        iMin = np.argmin(fseeFilt[N:500])
        if opts['idxQRmin'] == 'auto':
            idxMin = (round(N-N/2+iMin), round(N+N/2+iMin))
        else:
            idxMin = opts['idxQRmin'][iFile]
        
        # Select only interesting part of data
        iSel = fsee>0.5*np.max(fseeFilt)
        iStart = np.where(iSel)[0][0]
        iEnd = np.where(iSel)[0][-1]
        timeSel = time[iStart:iEnd]
        lmtcSel = lmtc[iStart:iEnd]
        fseeSel = fsee[iStart:iEnd]
        
        iEnd = np.where(np.gradient(fseeSel)>0)[0][-1]
        timeSel = timeSel[:iEnd]
        lmtcSel = lmtcSel[:iEnd]
        fseeSel = fseeSel[:iEnd]
        
        # Obtain the index at which the step occurs
        dlmtc = np.gradient(lmtcSel)
        iStep = np.argmin(dlmtc)
        
        # Obtain the indices slightly before and after the step
        iStop = np.argmax(np.gradient(fseeSel[iStep:]))+50
        # Obtain the index slightly before the step, this is based on:
            # As a consequence of the step lmtc quickly changes length, therefore 
            # we find the index just before the step as where dlmtc/dt is higher 
            # than 5% of the maximal value of dlmtc/dt during the step
        idxQRpreAvg = [idx for idx,var in enumerate(dlmtc[iStep-50:iStep]) if var>0.05*dlmtc[iStep]][-1]-50 
        # Obtain the index slightly after the step, this is based on: 
            # After the step is done, SEE force will increase, now we want to find
            # the last index in which SEE force is still decreasing
        idxQRpstAvg = [idx for idx,var in enumerate(np.gradient(fseeSel[iStep:iStep+iStop])) if var>0][0]
        # Obtain the indices slightly before and after the step  
        nPst = nQRsamp//2
        nPre = nQRsamp-nPst
        if opts['idxQRpre'] == 'auto':
            idxQRpre = (iStart+iStep+idxQRpreAvg-nPst, iStart+iStep+idxQRpreAvg+nPre)
        else:
            idxQRpre = opts['idxQRpre'][iFile]
        if opts['idxQRpst'] == 'auto':
            idxQRpst = (iStart+iStep+idxQRpstAvg-nPst, iStart+iStep+idxQRpstAvg+nPre)
        else:
            idxQRpst = opts['idxQRpst'][iFile]
        
        # Check it!
        idxs = {'fseeQRmin': idxMin, 'fseeQRpre': idxQRpre, 'fseeQRpst': idxQRpst}
        if 'dispFig' not in opts or opts['dispFig'] == True:
            idxs = DataPlot(time,lmtc,fsee,idxs,filename)
        
        dataQR[iFile]['time'] = time
        dataQR[iFile]['lmtc'] = lmtc
        dataQR[iFile]['fsee'] = fsee
        dataQR[iFile]['filename'] = filename
        dataQR[iFile]['idxQRmin'] = idxs['fseeQRmin']
        dataQR[iFile]['idxQRpre'] = idxs['fseeQRpre']
        dataQR[iFile]['idxQRpst'] = idxs['fseeQRpst']
    
    idxQRmin = [x['idxQRmin'] for x in dataQR]
    idxQRpre = [x['idxQRpre'] for x in dataQR]
    idxQRpst = [x['idxQRpst'] for x in dataQR]
    return dataQR, idxQRmin, idxQRpre, idxQRpst

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_qr_opts(opts):
    # Check user inputs
    if 'nQRsamp' in opts:
        if opts['nQRsamp'] != 'auto':
            if not isinstance(opts['nQRsamp'],int):
                raise TypeError(r"Invalid input. nQRsamp should be set to 'auto', or an integer number.")
            else:
                opts['nQRsamp'] = opts['nQRsamp']
        else:
            opts['nQRsamp'] = 5
    else:
        opts['nQRsamp'] = 5
    
    # Open files
    if 'dataDir' in opts:
        if not isinstance(opts['dataDir'],str):
            raise TypeError('Invalid input. dataDir should be a string.')
        else:
            dataDir = opts['dataDir']
            dataDir = dataDir + (not dataDir.endswith("/"))*'/'
            fQR = glob.glob(dataDir+r'*.csv')
            if len(fQR) < 1:
                raise ValueError("No data files found in data directory.")
    else:
        fQR = OpenFileDialog()
        
    nFiles = len(fQR)
    opts['nFiles'] = nFiles
    # Find how many columns are in te files
    nCol = 0
    for iFile,filepath in enumerate(fQR):
        df = pd.read_csv(filepath)
        if len(df.axes[1]) > nCol:
            nCol = len(df.axes[1])
    
    # Check user inputs
    if 'idxQRmin' in opts:
        if opts['idxQRmin'] != 'auto':
            if len(opts['idxQRmin']) != nFiles:
                raise ValueError('Invalid input. Your maually selected indices for fseeQRmin, but your indices do not match with the amount of data files. \n Number of data files = '+str(nFiles)+', and number of indices = '+str(len(opts['idxMin'])))
    else:
        opts['idxQRmin'] = 'auto'
    if 'idxQRpre' in opts:
        if opts['idxQRpre'] != 'auto':
            if len(opts['idxQRpre']) != nFiles:
                raise ValueError('Invalid input. Your maually selected indices for fseeQRpre, but your indices do not match with the amount of data files. \n Number of data files = '+str(nFiles)+', and number of indices = '+str(len(opts['idxMin'])))
    else:
        opts['idxQRpre'] = 'auto'
    if 'idxQRpst' in opts:
        if opts['idxQRpst'] != 'auto':
            if len(opts['idxQRpst']) != nFiles:
                raise ValueError('Invalid input. Your maually selected indices for fseeQRpst, but your indices do not match with the amount of data files. \n Number of data files = '+str(nFiles)+', and number of indices = '+str(len(opts['idxMin'])))
    else:
        opts['idxQRpst'] = 'auto'
    if 'iCols' in opts:
        if not isinstance(opts['iCols'], tuple):
            raise TypeError('Number of columns should be in a tuple.')
        elif not all(isinstance(x, int) for x in opts['iCols']):
            raise TypeError('Elements in lists should be integers.')
        elif len(opts['iCols']) != 3:
            raise ValueError('Number of colmuns should be 3 integers (time, lmtc and fsee).')
        else:
            opts['iCols'] = opts['iCols']
    else:
        opts['iCols'] = SelDataCol(3,nCol)
    return fQR,opts

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def sr(opts={},fSR={}):
    if fSR == {}:
        fSR,opts = get_sr_opts(opts)
    
    if len(fSR)<3:
        raise ValueError("Less than two datafiles are selected, at least 3 datafiles are necessary.") 
        return None, None, None, None
    
    iCols = opts['iCols']
    nFiles = opts['nFiles']
    nSRsamp = opts['nSRsamp']
    
    # Make listed-dic to store data   
    dataSR = [{'time': None, 'lmtc': None, 'fsee': None, 'filename': None, 'idxSRcon': None} for _ in range(nFiles)]
    # Now select data etc
    for iFile,filepath in enumerate(fSR):
        filename = filepath.rsplit('/', 1)[-1][:-4]
        data = pd.read_csv(filepath).T.to_numpy()
        time = data[iCols[0]]
        lmtc = data[iCols[1]] 
        fsee = data[iCols[2]]
        
        # First moving average filter (50 ms)
        dt = np.diff(time).mean()
        N = int(np.round(1/(20*dt)))
        B = np.ones(N)/N
        fseeFilt = signal.filtfilt(B,1,fsee)
        
        N = 20
        B = np.ones(N)/N
        lmtcFilt = signal.filtfilt(B,1,lmtc)
        
        # Select only interesting part of data
        fmax = np.max(fseeFilt)
        fmin = np.min(fseeFilt)
        iSel = fsee>0.5*(fmax-fmin)+fmin
        iStart = np.where(iSel)[0][0]
        iSel = fsee>0.01*(fmax-fmin)+fmin
        iEnd =  np.where(iSel)[0][-1]
        timeSel = time[iStart:iEnd]
        lmtcSel = lmtc[iStart:iEnd]
        fseeSel = fsee[iStart:iEnd]
        lmtcFiltSel = lmtcFilt[iStart:iEnd]
        # Now obtain the index at which the step occurs
        dlmtc       = np.gradient(lmtcSel)
        iStep       = np.argmin(dlmtc)

        # Select only interesting part of data (Part 2)
        timeSel = timeSel[iStep:]
        lmtcSel = lmtcSel[iStep:]
        fseeSel = fseeSel[iStep:]
        lmtcFiltSel = lmtcFiltSel[iStep:]
        vmtcSel = np.gradient(lmtcFiltSel)
        try:
            iEnd = np.argwhere(vmtcSel>0)[0][0]
        except:
            iEnd = len(vmtcSel)
        iEnd = int(0.7*iEnd) # should be in the first 70% of the ramp..
        
        # Select only interesting part of data (Part 3)
        timeSel = timeSel[:iEnd]
        lmtcSel = lmtcSel[:iEnd]
        fseeSel = fseeSel[:iEnd]
        
        # Obtain the index at which see force is most constant
        if opts['idxSRcon'] == 'auto':
            # Select 20 samples (=10ms) with most plateau
            fseeSlope = np.empty(len(fseeSel)-nSRsamp)*np.nan
            for iSample in range(0,len(fseeSlope)):
                iSel = np.arange(iSample,iSample+nSRsamp)
                fseeSlope[iSample] = np.polyfit(timeSel[iSel],fseeSel[iSel],1)[0]
            iEnd = np.nanargmin(np.gradient(fseeSlope))
            fseeSlope[iEnd+1:] = np.nan
            iMin = np.nanargmin(abs(fseeSlope))+iStart+iStep
            idxSRcon = (iMin,iMin+nSRsamp)
            
            fseeDot = np.gradient(fseeSel)
            fseeSlope = np.empty(len(fseeSel)-nSRsamp)*np.nan
            fseeDotAbsAvg = np.empty(len(fseeSel)-nSRsamp)*np.nan
            for iSample in range(0,len(fseeSlope)):
                iSel = np.arange(iSample,iSample+nSRsamp)
                fseeSlope[iSample] = np.polyfit(timeSel[iSel],fseeSel[iSel],1)[0]
                fseeDotAbsAvg[iSample] = np.abs(fseeDot[iSel]).mean()
            # iEnd = np.nanargmin(fseeSlope[20:])+20
            # fseeDotAbsAvg[iEnd+5:] = np.nan
            iMin = np.nanargmin(fseeDotAbsAvg)+iStart+iStep
            idxSRcon = (iMin,iMin+nSRsamp)
            
        else:
            idxSRcon = opts['idxSRcon'][iFile]
        
        # Check it!
        idxs = {'fseeSRcon': idxSRcon}
        if 'dispFig' not in opts or opts['dispFig'] == True:
            idxs = DataPlot(time,lmtc,fsee,idxs,filename)
          
        dataSR[iFile]['time'] = time
        dataSR[iFile]['lmtc'] = lmtc
        dataSR[iFile]['fsee'] = fsee
        dataSR[iFile]['filename'] = filename
        dataSR[iFile]['idxSRcon'] = idxs['fseeSRcon']
            
    idxSRcon = [x['idxSRcon'] for x in dataSR]
    return dataSR, idxSRcon

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_sr_opts(opts={}):
    # Check user inputs
    if 'nSRsamp' in opts:
        if opts['nSRsamp'] != 'auto':
            if not isinstance(opts['nSRsamp'],int):
                raise TypeError(r"Invalid input. nSRsamp should be set to 'auto', or an integer number.")
            else:
                opts['nSRsamp'] = opts['nSRsamp']
        else:
            opts['nSRsamp'] = 10
    else:
        opts['nSRsamp'] = 10
    
    # Open files
    if 'dataDir' in opts:
        if not isinstance(opts['dataDir'],str):
            raise TypeError('Invalid input. dataDir should be a string.')
        else:
            dataDir = opts['dataDir']
            dataDir = os.path.join(dataDir,'')
            fSR = glob.glob(dataDir+r'*.csv')
            if len(fSR) < 1:
                raise ValueError("No data files found in data directory.")
    else:
        fSR = OpenFileDialog()
    
    opts['nFiles'] = len(fSR)
    # Find how many columns are in te files
    nCol = 0
    for iFile,filepath in enumerate(fSR):
        df = pd.read_csv(filepath)
        if len(df.axes[1]) > nCol:
            nCol = len(df.axes[1])
    
    # Check user inputs
    if 'idxSRcon' in opts:
        if opts['idxSRcon'] != 'auto':
            if len(opts['idxSRcon']) != opts['nFiles']:
                raise ValueError('Invalid input. Your maually selected indices for fseeCon, but your indices do not match with the amount of data files. \n Number of data files = '+str(opts['nFiles'])+', and number of indices = '+str(len(opts['idxSRcon'])))
    else:
        opts['idxSRcon'] = 'auto'
    if 'iCols' in opts:
        if not isinstance(opts['iCols'], tuple):
            raise TypeError('Number of columns should be in a tuple.')
        elif not all(isinstance(x, int) for x in opts['iCols']):
            raise TypeError('Elements in lists should be integers.')
        elif len(opts['iCols']) != 3:
            raise ValueError('Number of colmuns should be 3 integers (time, lmtc and fsee).')
        else:
            opts['iCols'] = opts['iCols']
    else:
        opts['iCols'] = SelDataCol(3,nCol)
    
    return fSR,opts

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def isom(opts={}):
    # Check user inputs
    if 'durStimOffset' in opts:
        if opts['durStimOffset'] != 'auto':
            if not isinstance(opts['durStimOffset'],float):
                raise TypeError(r"Invalid input. durStimOffset should be set to 'auto', or a float.")
            else:
                durStimOffset = opts['durStimOffset']
        else:
            durStimOffset = 0.1
    else:
        durStimOffset = 0.1
    
    # Open files
    if 'dataDir' in opts:
        if not isinstance(opts['dataDir'],str):
            raise TypeError('Invalid input. dataDir should be a string.')
        else:
            dataDir = opts['dataDir']
            dataDir = os.path.join(dataDir,'')
            fISOM = glob.glob(dataDir+r'*.csv')
            if len(fISOM) < 1:
                raise ValueError("No data files found in data directory.")
    else:
        fISOM = OpenFileDialog()
    
    nFiles = len(fISOM)
    # Find how many columns are in te files
    nCol = 0
    for iFile,filepath in enumerate(fISOM):
        df = pd.read_csv(filepath)
        if len(df.axes[1]) > nCol:
            nCol = len(df.axes[1])
    
    # Check user inputs
    if 'iCols' in opts:
        if not isinstance(opts['iCols'], tuple):
            raise TypeError('Number of columns should be in a tuple.')
        elif not all(isinstance(x, int) for x in opts['iCols']):
            raise TypeError('Elements in lists should be integers.')
        elif len(opts['iCols']) != 4:
            raise ValueError('Number of colmuns should be 4 integers (time, lmtc, fsee and stim).')
        else:
            iCols = opts['iCols']
    else:
        iCols = SelDataCol(4,nCol)
        
    # Make listed-dic to store data   
    dataISOM = [{'time': None, 'lmtc': None, 'fsee': None, 'filename': None, 'idxSEL': None} for _ in range(nFiles)]
    # Now select data etc
    for iFile,filepath in enumerate(fISOM):
        filename = filepath.rsplit('/', 1)[-1][:-4]
        data = pd.read_csv(filepath).T.to_numpy()
        time = data[iCols[0]]
        lmtc = data[iCols[1]] 
        fsee = data[iCols[2]]
        stim = data[iCols[3]]
        
        _,tStimOn,tStimOff = get_stim(time,stim)
        if len(tStimOn)>1 or len(tStimOff)>1:
            print('Multiple stimulation trains detected, please manually select data interval.')
            idxSEL = (0,100)
        else:
            tStart = tStimOn
            tStop  = tStimOff+durStimOffset
            iStart = int(np.argmin(abs(time-tStart)))
            iStop  = int(np.argmin(abs(time-tStop)))
            idxSEL = (iStart,iStop)
        
        # Check it!
        idxs = {'fseeSEL': idxSEL}
        if 'dispFig' not in opts or opts['dispFig'] == True:
            idxs = DataPlot(time,lmtc,fsee,idxs,filename)
          
        dataISOM[iFile]['time'] = time
        dataISOM[iFile]['lmtc'] = lmtc
        dataISOM[iFile]['fsee'] = fsee
        dataISOM[iFile]['stim'] = stim
        dataISOM[iFile]['filename'] = filename
        dataISOM[iFile]['idxSEL']   = idxs['fseeSEL']
            
    idxSEL = [x['idxSEL'] for x in dataISOM]
    return dataISOM, idxSEL