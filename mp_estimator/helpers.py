# -*- coding: utf-8 -*-
"""
    helpers.py contains 'helper' functions, necessary for the parameter value
    estimation process. 

    author: Edwin D.H.M. Reuvers
    v1.0.0, september 2025
"""

#%% Import packages
import numpy as np
from scipy import optimize

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_stim(time,signal):
    """
    Finds the start and stop times of stimulation pulse trains in a signal.

    Inputs:
    - signal: 1D array-like stimulation over time (e.g., list or numpy array)

    Outputs:
    - pulse_trains: List of tuples where each tuple contains (start, stop) indices of a pulse train
    """
    
    # Determine an appropriate threshold as halfway between min and max of the signal
    threshold = (np.max(signal) + np.min(signal)) / 2
    
    # Convert signal to binary form based on the threshold
    binary_signal = np.where(signal > threshold, 1, 0)
    
    # Find the start and stop indices of all pulses
    iStartPulse = np.where(np.diff(binary_signal, prepend=0) == 1)[0]
    iStopPulse = np.where(np.diff(binary_signal, prepend=0) == -1)[0]
    
    if len(iStartPulse) == 0 or len(iStopPulse) == 0:
        return []
    
    # Estimate the minimum gap between pulse trains by finding the most common pulse width
    PulseWidths = iStopPulse - iStartPulse
    PulseWidthMed = np.median(PulseWidths)
    
    # Separate the start and stop indices into different pulse trains
    iStartTrain = [int(iStartPulse[0])]  # Convert to Python int
    iStopTrain = []
        
    for i in range(1, len(iStartPulse)):
        # If the gap between current stop and next start is too large, close the current train
        if iStartPulse[i] - iStopPulse[i-1] > PulseWidthMed:
            iStopTrain.append(int(iStopPulse[i-1]))  # Convert to Python int
            iStartTrain.append(int(iStartPulse[i]))  # Convert to Python int
    
    # Add the final pulse train stop
    iStopTrain.append(int(iStopPulse[-1]))  # Convert to Python int
        
    stimModel = get_block_signal(signal,iStartTrain, iStopTrain)
    tStimOn = time[iStartTrain]
    tStimOff = time[iStopTrain]
    return stimModel, tStimOn, tStimOff

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_block_signal(signal, pulse_train_starts, pulse_train_stops):
    """
    Creates a block-shaped signal around the pulse trains.
    
    Iputs:
    - signal: 1D array-like signal over time (e.g., list or numpy array)
    - pulse_train_starts: List of start indices of the pulse trains (Python integers)
    - pulse_train_stops: List of stop indices of the pulse trains (Python integers)
    
    Outputs:
    - block_signal: 1D numpy array of the same length as the input signal with block shapes (height 1) around pulse trains
    """
    
    # Initialize the block signal as a zero array with the same length as the input signal
    block_signal = np.zeros_like(signal)
    
    # Set the block signal to 1 between each pair of start and stop indices
    for start, stop in zip(pulse_train_starts, pulse_train_stops):
        block_signal[start:stop] = 1
    
    return block_signal

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def findModelFV(vceData, fceData, lcerelData, qData, muspar):
    """
    Estimates FCE and VCE on CE force-velocity relationship by minimising
    distance those vlaues and FCE and VCE data.
    
    Parameters:
        vceData : np.ndarray - measured VCE data
        fceData : np.ndarray - measured FCE data
        lcerelData : np.ndarray - normalized muscle length data
        qData : np.ndarray - active state data
        muspar : dict - muscle parameters (must contain 'a', 'b', 'fmax')
        
    Returns:
        vceEst : np.ndarray - predicted VCE
        fceEst : np.ndarray - predicted FCE
    """
    
    # Import packages
    from .hillmodel import Fce2Vce
    
    #
    n = len(vceData)
    vceEst = np.empty(n)
    fceEst = np.empty(n)
    
    # Precompute scaling factors
    scale_v = muspar['b'] / muspar['a'] * muspar['fmax']
    scale_f = muspar['fmax']
    
    for i in range(n):
        # Define distance function to minimize
        def distance(fce):
            vce = Fce2Vce(fce, qData[i], lcerelData[i], muspar)[0]
            vd = (vce - vceData[i]) / scale_v
            fd = (fce - fceData[i]) / scale_f
            return np.sqrt(vd**2 + fd**2)
        
        # Try fast 1D bounded optimizatio
        res = optimize.minimize_scalar(distance, bounds=(0, muspar['fmax']), method='bounded')
        
        # If this does not work: fallback to Nelder-Mead with initial guess
        if not res.success:
            # Formula initial guess first
            fcePossible = np.linspace(0, muspar['fmax'], 500)
            vcePossible = np.array([Fce2Vce(fce, qData[i], lcerelData[i], muspar)[0] for fce in fcePossible])
            vd = (vcePossible - vceData[i]) / scale_v
            fd = (fcePossible - fceData[i]) / scale_f
            d = np.sqrt(vd**2 + fd**2)
            fce0 = fcePossible[np.argmin(d)]  # initial guess
            
            res_nm = optimize.minimize(distance, fce0, method='Nelder-Mead', options={'xatol':1e-9})
            fceEst[i] = res_nm.x[0]  # res.x is array for Nelder-Mead
        else:
            fceEst[i] = res.x
        
        # Compute VCE once from optimized FCE
        vceEst[i] = Fce2Vce(res.x, qData[i], lcerelData[i], muspar)[0]
    
    return vceEst, fceEst

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def exd_guess(ACTdata, muspar):
    """
    Estimate an initial guess of the calcium dynamics activation and 
    deactivation time constant based on measured force
    
    Inputs:
        ACTdata: dict - containing data of used to estimate excitation dynamics parameters
        
    Outputs:
        tact : float - initial guess for calcium dynamics activation time constant
        tdeact : float - initial guess for calcium dynamics deactivation time constant
    """
    
    # Import packages
    from .objectives import cstfnc_act
    
    # Make dictionaries of only the part where there is stimulation
    ONdata = {}
    OFFdata = {}
    for iFile in ACTdata:
        time = ACTdata[iFile]['time']
        lmtc = ACTdata[iFile]['lmtc']
        fsee = ACTdata[iFile]['fseeData']
        _,t_off = ACTdata[iFile]['tStim']
        i_off = np.argmin(np.abs(time - t_off))
        
        ONdata[iFile] = {}
        ONdata[iFile]['time']       = time[0:i_off]
        ONdata[iFile]['lmtc']       = lmtc[0:i_off]
        ONdata[iFile]['fseeData']   = fsee[0:i_off]
        ONdata[iFile]['tStim']      = ACTdata[iFile]['tStim']
        ONdata[iFile]['gamma0']     = ACTdata[iFile]['gamma0']
        ONdata[iFile]['lcerel0']    = ACTdata[iFile]['lcerel0']

    # Find intial guess of calcium activation time-const first
    fun = lambda x: cstfnc_act(x,ONdata,muspar)[0]
    t0 = np.linspace(1e-3,1e-1,5)
    err = [fun([tOn,0.1]) for tOn in t0]
    tact = t0[np.argmin(err)]
    # Interpolate 
    coef = np.polyfit(t0,err,2)
    t0_int = np.linspace(t0[0],t0[-1],100)
    err_int = np.polyval(coef,t0_int)
    tact = t0_int[np.argmin(err_int)]
            
    # Then intial guess of calcium deactivation time-const first   
    fun = lambda x: cstfnc_act(x,ACTdata,muspar)[0]
    t0 = np.linspace(1e-3,1e-1,5)
    err = [fun([tact,tOn]) for tOn in t0]
    # Interpolate 
    coef = np.polyfit(t0,err,2)
    t0_int = np.linspace(t0[0],t0[-1],100)
    err_int = np.polyval(coef,t0_int)
    tdeact = t0_int[np.argmin(err_int)]
    
    return tact, tdeact
