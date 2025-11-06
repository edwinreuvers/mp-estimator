# -*- coding: utf-8 -*-
"""
    objectives.py contains all objective functions, that is: the objectives
    that are minised in order to estimate the MTC parameter values.

    author: Edwin D.H.M. Reuvers
    v1.0.0, september 2025
"""

#%% Import packages
import numpy as np

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def cstfnc_see(x,lseeDeltaStep,fseeQRpreData,fseeQRpstData,muspar):
    """
    cstfnc_see Computes rmse between data and modelled force of SEE 
    force-length relationship.
    
    Inputs:
        x                   =   array containing with on index 0 the stiffness  
                                shape constant, and index 1-end cSEE
        lseeDeltaStep       =   change in SEE length due to the quick-release
        fseeQRpreData       =   data SEE force before the step quick-release
        fseeQRpstData       =   data SEE force after the step quick-release
    
    Outputs:
        rmse                =   RMSE between data and modelled SEE force
        fseeQRpreModel      =   model SEE force before the quick-release
        fseeQRpstModel      =   model SEE force after the quick-release
    """
    
    #%% Read out muscle parameter values
    n       = muspar['n'] # exponential of SEE-force relationship
    sSEE    = x[0]  # [N/m^2] SEE stiffness shape constant
    cSEE    = x[1:] # [m]
    
    #%% Comute model SEE
    # Create new variables, such that model SEE force cannot be negative
    y = cSEE # [m]
    y[y<0] = 0 # [m]
    x = cSEE+lseeDeltaStep # [m]
    x[x<0] = 0 # [m]
    
    # Compute model SEE force before and after the quick-release
    fseeQRpreModel = sSEE*(y)**n # [N]
    fseeQRpstModel = sSEE*(x)**n # [N]
    
    #%% Compute rmse between data and modelled force
    e1 = (fseeQRpreData-fseeQRpreModel)**2
    e2 = (fseeQRpstData-fseeQRpstModel)**2
    rmse = (np.sum(e1) + np.sum(e2))**0.5
    
    #%% Output
    return rmse, fseeQRpreModel, fseeQRpstModel

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def cstfnc_pee(x,lpeePLUSlsee0,fpeeData,muspar):
    """
    cstfnc_pee Computes rmse between data and modelled force of PEE 
    force-length relationship.
    
    Inputs:
        x                   =   ..
        lpeeData            =   ..
        fpeeData            =   ..
    
    Outputs:
        rmse                =   RMSE between data and modelled PEE force
        fpeeModel           =   model PEE force
    """
    
    #%% Read-out muscle parameter values
    kpee    = x[0]  # shape parameter of EE curve
    cPEE    = x[1]  # 
    
    #%% Computations    
    x = lpeePLUSlsee0-cPEE
    x[x<0] = 0
    fpeeModel = kpee*(x)**2
    
    #%% Compute rmse between data and modelled force
    e = (fpeeData-fpeeModel)**2
    rmse = np.sum(e)**0.5
    
    #%% Outputs
    return rmse,fpeeModel

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def cstfnc_mtc(x,lmtcData,fseeData,parms,sa=0):
    """
    cstfnc_mtc Computes rmse between data and modelled force of multiple 
    trials
    
    Inputs:
        x                   =   array containing force-length parameter
                                    values (i.e., fmax, lce_opt, lsee0 and lpee0)
        data                =   dict containing the experimental data
                                    (e.g., lmtc[t], fsee[t] and stim[t])
        parms               =   dict containing the (muscle) parameter values
    
    Outputs:
        rmse                =   RMSE between data and modelled SEE force
    """
    
    #%% Import packages
    from .hillmodel import ForceEQ
    
    #%% Read-out muscle parameter values
    muspar = parms.copy()
    if sa == 1: # sensitivity for lsee0
        muspar['fmax']      = x[0];
        muspar['lce_opt']   = x[1]
        muspar['lpee0']     = muspar['cPEE']-muspar['lsee0']
    elif sa == 2: # sensitivity for lce_opt
        muspar['fmax']      = x[0];
        muspar['lsee0']     = x[1]
        muspar['lpee0']     = muspar['cPEE']-muspar['lsee0']
    elif sa == 3: # sensitivity for fmax
        muspar['lce_opt']   = x[0]
        muspar['lsee0']     = x[1]
        muspar['lpee0']     = muspar['cPEE']-muspar['lsee0']
    else:
        muspar['fmax']      = x[0];
        muspar['lce_opt']   = x[1]
        muspar['lsee0']     = x[2]
        muspar['lpee0']     = muspar['cPEE']-muspar['lsee0']
        
    #%% Compute model SEE force
    fseeModel = ForceEQ(lmtcData,1,muspar)[0]
    
    #%% Compute rmse between data and modelled force
    e = (fseeData-fseeModel)**2
    rmse = np.sum(e)**0.5
    
    #%% Outputs
    return rmse, fseeModel

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def cstfnc_fv(x,vceData,fceData,lcerelData,parms,sa=0):
    """
    cstfnc_fv Computes rmse between data and modelled force of multiple 
    trials
    
    Inputs:
        x                   =   array containing force-velocity parameter
                                    values (i.e., arel & brel)
        data                =   dict containing the experimental data
                                    (e.g., lmtc[t], fsee[t] and stim[t])
        parms               =   dict containing the (muscle) parameter values
    
    Outputs:
        err                =   measure of differece between data and modelled 
                                    SEE force and/or data and modelled CE
                                    velocity
    """
    
    #%%
    from .hillmodel import ActState
    from .helpers import findModelFV
    
    #%% Read-out muscle parameter values
    muspar = parms.copy()
    if len(x) == 2:
        muspar['a'] = x[0]
        muspar['b'] = x[1]
    elif len(x) == 3:
        muspar['a'] = x[0]
        muspar['b'] = x[1]
        muspar['fmax'] = x[2]
    else:
        raise TypeError("")
        
    #%% Computations
    qData = ActState(1,lcerelData,muspar)[0]  # [ ] active state
    vceModel,fceModel = findModelFV(vceData,fceData,lcerelData,qData,muspar)
    
    # Compute error
    d1  = (fceData-fceModel)
    d2  = (vceData-vceModel)
    s1  = np.abs(np.max(fceData)-np.min(fceData))
    s2  = np.abs(np.max(vceData)-np.min(vceData))
    e1  = (d1/s1)**2
    e2  = (d2/s2)**2
    err = np.sum(e1) + np.sum(e2) + sum(vceModel>0)*1        
    
    #%% Outputs
    return err,fceModel,vceModel

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def cstfnc_act(x,data,parms):
    """
    cstfnc_act Computes rmse between data and modelled force of multiple 
    trials
    
    Inputs:
        x                   =   array containing activation and deactivation
                                time constant
        data                =   dict containing the experimental data
                                    (e.g., lmtc[t], fsee[t] and stim[t])
        parms               =   dict containing the (muscle) parameter values
    
    Outputs:
        rmsd                =   RMSD between data and estimated SEE force
    """
    
    #%% Import packages  
    from .hillmodel import SolveSimuMTC
    def compute_rmsd(x,y):
        n = ((x-y)**2).mean()
        return n**(1/2)
    
    #%% Read-out muscle parameter values
    muspar = parms.copy()
    muspar['tact'] = x[0]
    muspar['tdeact'] = x[1]
    
    #%% Compuations    
    rmsd = []
    for iFile,dc in enumerate(data):
        _,solstr = SolveSimuMTC(data[dc]['gamma0'],data[dc]['lcerel0'],muspar,data[dc])
        fseeData = data[dc]['fseeData']
        fseeModel = solstr[9]
        
        rmsd.append(compute_rmsd(fseeData,fseeModel))
        data[dc]['fseeModel'] = fseeModel
    
    #%% Outputs
    return np.sum(rmsd), data