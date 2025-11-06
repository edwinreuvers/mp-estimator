# -*- coding: utf-8 -*-
"""
    estimate.py contains all functions to estimate the muscle parameter values

    author: Edwin D.H.M. Reuvers
    v1.0.0, september 2025
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from .helpers import exd_guess, get_stim
from .hillmodel import ActState, Fce2Vce, ForceEQ, ForceLength, GetForceCE, SimuMTC, SolveSimuMTC
from .objectives import cstfnc_act, cstfnc_fv, cstfnc_mtc, cstfnc_pee, cstfnc_see

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def im(dataQR,dataSR,parin,do_print=False):
    muspar = parin.copy()
    
    #%% Step 1: Estimate muspar for the firs time
    estpar = fl(dataQR,muspar)[0]
    estpar = fv(dataSR,muspar,estpar)[0]
    
    #%% Step 2: Iterative process
    iRound = 1
    maxDifference = 1
    while maxDifference > 0.001:
        oldpar = estpar.copy()
        maxDifference = 0
        estpar = fl(dataQR,muspar,estpar)[0]
        estpar = fv(dataSR,muspar,estpar)[0]
        for parName in oldpar:
            # Compute relative difference
            difference = np.max(np.abs((oldpar[parName] - estpar[parName]) / oldpar[parName]))

            # Check if it's the largest so far
            if difference > maxDifference:
                maxDifference = difference
                maxDiffKey = parName

            # Smooth update
            estpar[parName] = oldpar[parName] * 0.7 + 0.3 * estpar[parName]
        
        if do_print == True:
            print(f"For round no. {iRound}: max. difference = {maxDifference * 100:.4f}% in parameter '{maxDiffKey}'")
        iRound = iRound+1
        
    estpar,dataQRout = fl(dataQR,muspar,estpar)
    estpar,dataSRout = fv(dataSR,muspar,estpar)
    return estpar, dataQRout, dataSRout

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def fl(dataQR,muspar,estparr={}):
    estpar = estparr.copy()
    
    #%% Step 0: Get data variables    
    nFiles = len(dataQR)
    fseeQRmin, fseeQRpre, fseeQRpst, lmtcQRmin, lmtcQRpre, lmtcQRpst, lceDeltaStep = \
    [np.empty(nFiles) * np.nan for _ in range(7)]
    for iFile,data in enumerate(dataQR):
        time        = data['time']
        lmtc        = data['lmtc']
        fsee        = data['fsee']
        filename    = data['filename']
        
        idxQRmin = data['idxQRmin']
        idxQRpre = data['idxQRpre']
        idxQRpst = data['idxQRpst']
        
        # Now assign the variable which we use to estimate musparms!
        fseeQRmin[iFile] = fsee[idxQRmin[0]:idxQRmin[1]].mean()
        fseeQRpre[iFile] = fsee[idxQRpre[0]:idxQRpre[1]].mean()
        fseeQRpst[iFile] = fsee[idxQRpst[0]:idxQRpst[1]].mean()
        lmtcQRmin[iFile] = lmtc[idxQRmin[0]:idxQRmin[1]].mean()
        lmtcQRpre[iFile] = lmtc[idxQRpre[0]:idxQRpre[1]].mean()
        lmtcQRpst[iFile] = lmtc[idxQRpst[0]:idxQRpst[1]].mean()

        # Get CE shortening during the step
        if 'a' in estpar and 'b' in estpar: # we can only do this if FV parms are known!
            # Inputs of simulation
            strct = {}
            strct['time'] = time[idxQRpre[0]:idxQRpst[1]]-time[idxQRpre[0]]       
            strct['lmtc'] = lmtc[idxQRpre[0]:idxQRpst[1]]
            strct['tStim'] = [strct['time'][0], strct['time'][-1]]
            
            # Values at t=0
            lmtc0 = strct['lmtc'][0]
            gamma0 = 1
            lcerel0 = ForceEQ(lmtc0,gamma0,{**muspar,**estpar})[1]
            state0 = [gamma0, lcerel0]
            
            # Do simulation
            tspan = strct['tStim']
            fun = lambda t, x: SimuMTC(t,x,{**muspar,**estpar},strct)[0:2]
            sol = integrate.solve_ivp(fun,tspan,state0,method='RK45',max_step=5e-4,rtol=1e-4,atol=1e-8,t_eval=strct['time'])
            
            # Compute CE length change as a result of the step
            lce = sol.y[1]*estpar['lce_opt'] # [m]
            lceQRpre = lce[:np.diff(idxQRpre)[0]]
            lceQRpst = lce[-np.diff(idxQRpst)[0]:]
            lceDeltaStep[iFile] = lceQRpst.mean() - lceQRpre.mean()
        else: # During the first run, we do not know the FV parms yet, so we cannot compute CE length change
            lceDeltaStep[iFile] = 0
    
    lmtcDeltaStep = lmtcQRpst-lmtcQRpre # [m] change in lmtc due to quick-release
    lseeDeltaStep = lmtcDeltaStep - lceDeltaStep  # [m] change in lsee due to quick-release
    
    #%% Step 1: Estimate SEE stifness parameter (ksee)
    # Compute value of ksee for each individual QR
    kseeIndv = ((np.sqrt(fseeQRpst) - np.sqrt(fseeQRpre))/(lseeDeltaStep))**2# [N/m.^2]
    # The average is a good estimate of the temporary SEE stiffness parameter ksee
    ksee0 = np.mean(kseeIndv) # [N/m.^2] This is then our initial guess for the optimisation
    cSEE0 = (fseeQRpre/ksee0)**(1/2) # [m]
    seePar0 = np.insert(cSEE0,0,ksee0)
    fun = lambda x: cstfnc_see(x,lseeDeltaStep,fseeQRpre,fseeQRpst,muspar)[0]
    result = optimize.minimize(fun,seePar0,method='Nelder-Mead',options={'maxiter':1e7})
    if result.success !=True:
        seePar = result.x
        _, fseeQRpreModel, fseeQRpstModel = cstfnc_see(seePar,lseeDeltaStep,fseeQRpre,fseeQRpst,muspar)
        plt.plot(seePar[1:],fseeQRpre,'.',seePar[1:],fseeQRpreModel,'.')
        plt.plot(seePar[1:]+lseeDeltaStep,fseeQRpst,'.',seePar[1:]+lseeDeltaStep,fseeQRpstModel,'.')
        breakpoint()
    seePar = result.x.tolist()
    estpar['ksee'] = seePar[0]
    
    #%% Step 2: Estimate PEE curve
    esee = (((fseeQRmin>0)*fseeQRmin+(fseeQRmin<=0)*0)/estpar['ksee'])**(1/2)
    if 'lsee0' and 'fmax'  and 'lce_opt' in estpar:
        lce = lmtcQRmin-esee-estpar['lsee0']
        lcerel = lce/estpar['lce_opt']
        fceQRmin = ForceLength(lcerel,muspar)[0]*muspar['q0']*estpar['fmax'] 
    else:
        fceQRmin = 0
    fpeeQRmin = fseeQRmin-fceQRmin 
    fpeeQRmin = (fpeeQRmin>0)*fpeeQRmin + (fpeeQRmin<0)*0 # [N] PEE force cannot be smaller than 0!    
    lpeePLUSlsee0 = lmtcQRmin - (((fseeQRmin>0)*fseeQRmin+(fseeQRmin<=0)*0)/estpar['ksee'])**(1/2) # [m] PEE length + SEE slack length
    coef = np.polyfit(lpeePLUSlsee0,fseeQRmin,2)
    coefd = np.polyder(coef,1)
    kpeeGuess = coefd[0]/2
    lpee0Guess = -coefd[1]/coefd[0]
       
    peePar0 = np.array([kpeeGuess,lpee0Guess])
    fun = lambda x: cstfnc_pee(x,lpeePLUSlsee0,fpeeQRmin,{**muspar,**estpar})[0]
    result = optimize.minimize(fun,peePar0,method='Nelder-Mead')
    if result.success !=True:
        fpeeQRminMdl = cstfnc_pee(peePar0,lmtcQRmin,fseeQRmin,{**muspar,**estpar})[1]
        plt.plot(lmtcQRpre - (fseeQRmin/estpar['ksee'])**(1/2),fseeQRmin,'.',lmtcQRpre - (fseeQRmin/estpar['ksee'])**(1/2),fpeeQRminMdl,'.')
        breakpoint()
    peePar = result.x.tolist()
    estpar['kpee'],estpar['cPEE'] = peePar
        
    #%% Step 3: Estimate MTC curve  
    fmax0 = fseeQRpre.max()
    iMaxFile = np.argmax(fseeQRpre)
    lcerel0 = (np.sqrt(muspar['w']**2-fseeQRpre/fmax0*muspar['w']**2))*np.sign(lmtcQRpre-lmtcQRpre[iMaxFile])+1
    eseeBeforeStep = np.sqrt(fseeQRpre/estpar['ksee'])
    A = np.array([lcerel0, np.ones(np.size(lcerel0))])
    b = lmtcQRpre-eseeBeforeStep
    x0 = np.linalg.lstsq(A.T,b.T,rcond=None)[0] # Solve for lsee0+lcerel*lceopt = lmtc-esee
    if result.success !=True:
        breakpoint()
    cePar0 = np.array([fmax0, x0[0], x0[1]])
    fun = lambda x: cstfnc_mtc(x,lmtcQRpre,fseeQRpre,{**muspar,**estpar})[0]
    result = optimize.minimize(fun,cePar0,method='Nelder-Mead')
    cePar = result.x.tolist()
    estpar['fmax'],estpar['lce_opt'],estpar['lsee0'] = cePar

    #%% Step 4: Calculate muscle parameters from temporary 'parameters'
    estpar['lpee0'] = estpar['cPEE']-estpar['lsee0']
    del estpar['cPEE']
    
    #%% Step 5: Compute output variables
    lseeQRpre    = np.array(seePar[1:])+estpar['lsee0']
    lseeQRpst    = lseeQRpre+lseeDeltaStep
    lpeeQR       = lmtcQRpre - (((fseeQRmin>0)*fseeQRmin+(fseeQRmin<=0)*0)/estpar['ksee'])**(1/2) - estpar['lsee0']
    fpeeQR       = fseeQRmin
    
    dataOut = {}
    dataOut['lseeQRpre'] = lseeQRpre
    dataOut['lseeQRpst'] = lseeQRpst
    dataOut['fseeQRpre'] = fseeQRpre
    dataOut['fseeQRpst'] = fseeQRpst
    dataOut['lpeeQR'] = lpeeQR
    dataOut['fpeeQR'] = fpeeQR
    dataOut['lmtcQRpre'] = lmtcQRpre
        
    return estpar, dataOut

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def fv(dataSR,muspar,estparr={}):
    estpar = estparr.copy()
    
    #%% Step 0: Get data variables    
    nFiles = len(dataSR)
    vce, fce, lcerel = [np.empty(nFiles) * np.nan for _ in range(3)]
    for iFile,data in enumerate(dataSR):
        time        = data['time']
        lmtc        = data['lmtc']
        fsee        = data['fsee']
        filename    = data['filename']
        
        idxSRcon = data['idxSRcon']
        
        # Now assign variables during plateau of ramp
        iRamp = slice(idxSRcon[0], idxSRcon[1])
        timeRamp = time[iRamp]
        lmtcRamp = lmtc[iRamp]
        fseeRamp = fsee[iRamp]

        # Description..
        if 'ksee' in estpar:
            fceRamp,_,_,lceRamp = GetForceCE(lmtcRamp,fseeRamp,{**muspar,**estpar})[0:4]
            lcerelRamp = lceRamp/estpar['lce_opt']
        else:
            fceRamp     = fseeRamp
            lceRamp     = lmtcRamp
            lcerelRamp  = lmtcRamp/lmtcRamp
        
        # Store for computations
        vce[iFile]      = np.polyfit(timeRamp,lceRamp,1)[0]
        fce[iFile]      = np.mean(fceRamp)
        lcerel[iFile]   = np.mean(lcerelRamp)
    
    #%% Step 1: Compute FV variables
    q           = ActState(1,lcerel,muspar)[0]  # [ ] active state
    fisomrel    = ForceLength(lcerel,muspar)[0] # [ ] isometric CE force
    
    if 'fmax' in estpar:
        fmax = estpar['fmax']
        # Find initial guess of arel &  brel, by finding the lest-squares solution to:
        # - a*q*vce + b*(fce-q*fisomrel*fmax) = fce*vce    
        A = np.array([-q*vce, fce-q*fisomrel*fmax]).T
        b = fce*vce
        fvPar0 = np.linalg.lstsq(A, b, rcond=None)[0]
        # fvPar0 = np.insert(fvPar0,2,fmax)
        fun = lambda x: cstfnc_fv(x,vce,fce,lcerel,{**muspar,**estpar})[0]
        result = optimize.minimize(fun,fvPar0,method='Nelder-Mead',bounds=((1e-6, None),(1e-6, None)),options={'disp': False,'maxiter':1e4})
        if result.success !=True:
            breakpoint()
        fvPar = result.x.tolist()
        estpar['a'],estpar['b'] = fvPar
        _,fceModel,vceModel = cstfnc_fv(fvPar,vce,fce,lcerel,{**muspar,**estpar})
    else:
        fmax0 = fsee.max()
        # Find initial guess of arel &  brel, by finding the lest-squares solution to:
        # - a*q*vce + b*(fce-q*fisomrel*fmax) = fce*vce    
        A = np.array([-q*vce, fce-q*fisomrel*fmax0]).T
        b = fce*vce
        fvPar0 = np.insert(np.linalg.lstsq(A, b, rcond=None)[0],2,fmax0)
        # fvPar0 = np.insert(fvPar0,2,fmax)
        fun = lambda x: cstfnc_fv(x,vce,fce,lcerel,{**muspar,**estpar})[0]
        result = optimize.minimize(fun,fvPar0,method='Nelder-Mead',bounds=((1e-6, None),(1e-6, None),(1e-6, None)),options={'disp': False,'maxiter':1e4})
        if result.success !=True:
            breakpoint()
        fvPar = result.x.tolist()
        estpar['a'],estpar['b'],estpar['fmax'] = fvPar
        _,fceModel,vceModel = cstfnc_fv(fvPar,vce,fce,lcerel,{**muspar,**estpar})
    
    #%% Step 2: Compute output variables
    dataOut = {}
    dataOut['vceSR'] = vce
    dataOut['fceSR'] = fce
    dataOut['lcerelSR'] = lcerel
    
    return estpar, dataOut

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def act(dataISOM,defpar,estparr={},do_print=False):
    estpar = estparr.copy()
    
    #%% Step 0: Get data variables    
    ACTdata = {}
    for iFile,data in enumerate(dataISOM):
        time        = data['time']
        lmtc        = data['lmtc']
        fsee        = data['fsee']
        stim        = data['stim']
        filename    = data['filename']
        
        _,tStimOn, tStimOff = get_stim(time,stim)
        
        # Now assign variables during plateau of ramp
        idxSEL = data['idxSEL']
        idxSel = slice(idxSEL[0], idxSEL[1])
        tStimOn = tStimOn-time[idxSEL[0]]
        tStimOff = tStimOff-time[idxSEL[0]]
        time = time[idxSel]-time[idxSEL[0]]
        lmtc = lmtc[idxSel]
        fsee = fsee[idxSel]
        
        # Get initial state
        gamma0 = defpar['gamma_0']
        lcerel0 = ForceEQ(lmtc[0],gamma0,{**defpar,**estpar})[1]
        
        # Store in dict
        ACTdata[iFile] = {}
        ACTdata[iFile]['time'] = time
        ACTdata[iFile]['fseeData'] = fsee
        ACTdata[iFile]['lmtc'] = lmtc
        ACTdata[iFile]['tStim'] = [float(tStimOn[0]), float(tStimOff[0])]
        ACTdata[iFile]['gamma0'] = gamma0
        ACTdata[iFile]['lcerel0'] = lcerel0
    
    #%% Callback
    if do_print == True:
        class logger:
            def __init__(self):
                self.iteration = 0
    
            def callback(self, x):
                rmsd = fun(x)
                print(f"For iteration {self.iteration}: rmsd = {rmsd}")
                self.iteration += 1
        callback = logger().callback
    else:
        callback = None
    
    #%% Step 1: find initial guess
    if do_print == True:
        print('Trying to find initial guess for calcium dynamics time constants.')
    actPar0 = exd_guess(ACTdata,{**defpar,**estpar})
    if do_print == True:
        print(f'Initial guesses found. Tact = {actPar0[0]*1e3:0.1f} ms & tdeact = {actPar0[1]*1e3:0.1f} ms')
    #%% Step 2: Get activation dynamics parameters
    bnds = ((1e-3, 1), (1e-3, 1))
    fun = lambda x: cstfnc_act(x,ACTdata,{**defpar,**estpar})[0]   
    result = optimize.minimize(fun,actPar0,method='Nelder-Mead', options={'xatol': 1e-4}, callback=callback)
    actPar = result.x.tolist()
    estpar['tact'],estpar['tdeact'] = actPar
    _,ACTdata = cstfnc_act(actPar,ACTdata,{**defpar,**estpar})
    
    #%% Step 2: Compute output variables
    dataOut = [{k: v for k, v in ACTdata[d].items() if k not in {'gamma0', 'lcerel0'}} for d in ACTdata]   
    return estpar,dataOut
