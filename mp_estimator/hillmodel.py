# -*- coding: utf-8 -*-
"""
    hillmodel.py contains all functions regarding the Hill-type MTC model.

    author: Edwin D.H.M. Reuvers
    v1.0.0, september 2025
"""

#%% Import packages
import numpy as np
from scipy.integrate import solve_ivp, trapezoid

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def SimuMTC(t,state,parms,sol):
    """
    SimuMTC Computes state-derivatives (gammad and vcerel) for an Hill-type 
        MTC model.
    
    Inputs:
        t           =   time [s]
        state       =   state of the muscle model (gamma: the normalised
                            concentration Ca2+ between the filament and 
                            lcerel: relative CE length, i.e. CE length divided 
                            by optimum CE length) [ , ]
        parms       =   dict with muscle parameter values
        sol         =   dict with time-axis, lmtc(t) and either stim(t) or
                            stim onset- and offset times
    
    Outputs:
        gammad      =   time-derivative of gamma [1/s]
        vcerel      =   time-dertivative of relative CE length [1/s]
        y           =   list with other variables (i.e., lmtc, stim, q, etc.)
    """
    
    #%% Unravel parms & state
    gamma = state[0] # [ ]
    lcerel = state[1] # [ ]
    
    lmtc = np.interp(t,sol['time'],sol['lmtc']) # [m]
    
    #%% Get stim
    try: # first try if stim onset- and offset times are present
        tStim = sol['tStim'] # [s]
        tStim = np.atleast_2d(tStim)
        stim = np.zeros_like(lcerel)
        for start, end in tStim:
            stim[(t >= start) & (t <= end)] = 1
    except KeyError: # if stim onset- and offset time are not present stim(t) should be present
        stim = np.interp(t,sol['time'],sol['stim']) # [ ]
    
    #%% Activation dynamics
    gamma_0 = parms['gamma_0'] # [ ]
    gamma = (gamma>gamma_0)*gamma + (gamma<=gamma_0)*gamma_0 # [ ]
    q = ActState(gamma,lcerel,parms)[0] # [ ]
    gammad = (stim>=gamma)*((stim*(1-gamma_0)-gamma + gamma_0)/parms['tact']) + (stim<gamma)*((stim*(1-gamma_0)-gamma + gamma_0)/parms['tdeact']) # [1/s]
    
    #%% Get SEE, PEE and CE forces
    lce = lcerel*parms['lce_opt'] # [m]
    lsee = lmtc-lce # [m]
    lpee = lce # [m]
    fisomrel = ForceLength(lcerel, parms)[0] # [ ]
    fsee,fpee,_,_ = LEE2Force(lsee,lpee,parms) # [N, N, N/m, N/m]
    fce = fsee-fpee # [N]
    fcerel = fce/parms['fmax'] # [ ]
    
    #%% Get CE velocity  
    vcerel = Fce2Vce(fce,q,lcerel,parms)[1] # [1/s]
    
    #%% Debugging   
    if np.isnan(vcerel).any():
        breakpoint()
    
    #%% Output
    y = [lmtc, stim, q, lsee, lpee, fisomrel, fsee, fpee, fce, fcerel]
    
    return gammad, vcerel, y

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def SolveSimuMTC(gamma0,lcerel0,parms,strct):
    """
    SolveSimuMTC Perform a forward simulation for an Hill-type MTC model.
    
    Inputs:
        gamma0      =   the normalised concentration Ca2+ between the filament
                            at t=0 [ ]
        lcerel0     =   the relative CE length, i.e. CE length divided by
                            optimum CE length at t=0 [ ]
        parms       =   dict with muscle parameter values
        sol         =   dict with time-axis, lmtc(t) and either stim(t) or
                            stim onset- and offset times
    
    Outputs:
        Wmech       =   SEE mechanical work [J]
        y           =   list with variables (i.e., time, lmtc, stim, states, 
                            forces, etc.)
    """
        
    #%%
    state0 = [gamma0, lcerel0]
    tspan = [0, strct['time'][-1]]
    fun = lambda t, x: SimuMTC(t,x,parms,strct)[0:2]
    sol = solve_ivp(fun,tspan,state0,method='Radau',max_step=1e-3,rtol=1e-6,atol=1e-9,t_eval=strct['time'])
    time = sol.t; state = sol.y;
    gammad,vcerel,sol = SimuMTC(time,state,parms,strct)
    
    #%% Unravel the solution
    lmtc = sol[0]
    stim = sol[1]
    gamma = state[0]
    lcerel = state[1]
    q = sol[2]
    lsee = sol[3]
    lpee = sol[4]
    fisomrel = sol[5]
    fsee = sol[6]
    fpee = sol[7]
    fce = sol[8]
    fcerel = sol[9]
    
    y = [time, lmtc, stim, gamma, lcerel, q, lsee, lpee, fisomrel, fsee, fpee, fce, fcerel]
    
    #%% Output
    Wmech = -trapezoid(fsee,lmtc)
    return Wmech, y

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def ActState(gamma,lcerel,muspar):
    """
    ActState Computes the active state based on the relative amount of Ca2+ 
    between the filaments and the relative CE length. This version is based on 
    Hatze (1981), p. 37-41, but slightly modified such that the parameter
    values has physiological relevance and that it is easier to use in 
    Optimal Control.
    
    Inputs:
        gamma       =   the elative amount of Ca2+ between the filaments [ ]
        lcerel      =   the relative CE length (i.e. CE length divided by 
                            optimum CE length [ ]
        parms       =   dict with activation dynamics parameters (constants)
    
    Outputs:
        q           =   the active state (i.e. the relative amount of Ca2+ 
                        bound to troponin C) [ ]
        dqdlcerel   =   the partial dertivative of q with respect to lcerel [ ]
        gamma05 	=   the value of gamma where q=0.5 [ ]
    """
        
    #%% Unravel parameter values
    q0 = muspar['q0']; kCa = muspar['kCa']; a_act = muspar['a_act']
    b_act = muspar['b_act']; 
    
    #%% Computations
    a1_act = np.log10(np.exp(a_act))
    B_act = b_act[0]+b_act[1]*lcerel+b_act[2]*lcerel**2
    
    q = q0 + (1-q0) / (1 + (kCa*gamma)**(a1_act)*np.exp(a_act*B_act))
    dqdlcerel = (a_act*np.exp(a_act*(B_act))*(gamma*kCa)**(np.log(np.exp(a_act))/np.log(10))*(q0-1)*(b_act[1] + 2*b_act[2]*lcerel))/(np.exp(a_act*(B_act))*(gamma*kCa)**(np.log(np.exp(a_act))/np.log(10)) + 1)**2;
    gamma05 = np.where((lcerel <= 1 + muspar['w']) & (lcerel >= 1 - muspar['w']), ((1-0.5)/(kCa**a1_act*np.exp(a_act*B_act)*(0.5-q0)))**(1/a1_act), np.nan)
        
    #%% Output
    return q,dqdlcerel,gamma05

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def Fce2Vce(fce,q,lcerel,muspar):
    """
    Fce2Vce Computes the CE velocity based on the force-velocity relationship.
    
    Inputs:
        fce         =   the CE force  [N]
        q           =   the active state (i.e. the relative amount of Ca2+ 
                        bound to troponin C) [ ]
        lcerel      =   the relative CE length (i.e. CE length divided by 
                            optimum CE length [ ]
        muspar      =   dict with muscle parameter values
        
    Ouputs:
        vce         =   the CE velocity [m/s]
        vcerel      =   the relative CE velocity (i.e. CE velocity divided by 
                            optimum CE length) [1/s], will be only computed 
                            if optimum CE length is known, else output is NaN

    """
    
    #%% Unravel parameter values    
    a_c = muspar['a']; b_c = muspar['b']; fasymp = muspar['fasymp']; fmax = muspar['fmax']
    slopfac = muspar['slopfac']; vfactmin = muspar['vfactmin']; q0 = muspar['q0']  
    sloplin = vfactmin*b_c/(slopfac*0.005*0.0975*(fmax+a_c))
    
    #%% Scale arel and brel if necessary
    fisomrel = ForceLength(lcerel,muspar)[0]
    
    # Smooth version of KvS brel(q)
    q0_b = (np.log(1/vfactmin-1)+q0*22)/22
    b = b_c/(1+np.exp(-22*(q-q0_b)))
    
    a = a_c
    a = (lcerel>1)*a*fisomrel + (lcerel<=1)*a

    #%%
    dvdf_isom_con = b/(q*(fisomrel*fmax+a)) # slope in the isometric point at wrt concentric part
    dvdf_isom_ecc = dvdf_isom_con/slopfac # slope in the isometric point at wrt eccentric part
    dFdvcon0      = 1/dvdf_isom_con
    s_as          = 1/sloplin
    p1 = -(fisomrel*q*(fasymp*fmax - fmax))/(s_as - dFdvcon0*slopfac) 
    p2 =  (fisomrel**2*q**2*(fasymp*fmax - fmax)**2)/(s_as - dFdvcon0*slopfac)
    p3 =  -fasymp*fisomrel*q*fmax;
    p4 =  -s_as

    #%% Regions
    r_c1 = (((fce/q) <= fisomrel*fmax) & (dvdf_isom_con<=sloplin)) # Concentric, dvdf_isom_con<=sloplin (normal)	
    r_c2 = (((fce/q) <= fisomrel*fmax) & (dvdf_isom_con>sloplin)) # Concentric dvdf_isom_con>sloplin (defective case)
    r_e1 = (((fce/q) > fisomrel*fmax) & (dvdf_isom_ecc<=(sloplin/slopfac))) # Eccentric, dvdf_isom_ecc<=sloplin (normal) 
    r_e2 = (((fce/q) > fisomrel*fmax) & (dvdf_isom_ecc>(sloplin/slopfac))) # Eccentric, dvdf_isom_ecc>sloplin (defective case)
    
    #%% Compute CE velocity
    vce_c1 = (b*(fce-q*fisomrel*fmax)/(fce+q*a)) # Concentric, dvdf_isom_con<=sloplin (normal)	
    vce_c2 = (sloplin*(fce-q*fisomrel*fmax)) # Concentric dvdf_isom_con>sloplin (defective case)
    vce_e1 = ((-(fce + p3 + p1*p4 + (fce**2 - 2*fce*p1*p4 + 2*fce*p3 + p1**2*p4**2 - 2*p1*p3*p4 + p3**2 + 4*p2*p4)**(1/2))/(2*p4))) # Eccentric, dvdf_isom_ecc<=sloplin (normal) 
    vce_e2 = ((sloplin/slopfac)*(fce-q*fisomrel*fmax)) # Eccentric, dvdf_isom_ecc>sloplin (defective case)
    
    #%% Output
    vce = r_c1*vce_c1 + r_c2*vce_c2 + r_e1*vce_e1 + r_e2*vce_e2
    if 'lce_opt' in muspar:
        vcerel = vce/muspar['lce_opt']
    else:
        vcerel = None
    
    return vce, vcerel, [vce_c1, vce_c2, vce_e1, vce_e2], [r_c1, r_c2, r_e1, r_e2]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def Vce2Fce(vce,q,lcerel,muspar):
    """
    Vce2Fce Computes the CE force based on the force-velocity relationship.
    
    Inputs:
        vcerel      =   the relative CE velocity (i.e. CE velocity divided by 
                            optimum CE length [1/s]
        q           =   the active state (i.e. the relative amount of Ca2+ 
                            bound to troponin C) [ ]
        lcerel      =   the relative CE length (i.e. CE length divided by 
                            optimum CE length [ ]
        muspar      =   dict with muscle parameter values
        
    Ouputs:
        fcerel      =   the relative CE force (i.e., CE force normalised by
                            maximal isometric CE force) [ ]
        fce         =   the CE force [N]
    """
            
    #%% Unravel parameter values
    a_c = muspar['a']
    b_c = muspar['b']
    fasymp = muspar['fasymp']
    slopfac = muspar['slopfac']
    vfactmin = muspar['vfactmin']
    q0 = muspar['q0']
    fmax = muspar['fmax']
    
    sloplin = vfactmin*b_c/(slopfac*0.005*0.0975*(fmax+a_c))
    
    #%% Scale arel and brel if necessary
    fisomrel = ForceLength(lcerel,muspar)[0]
    
    # Smooth version of KvS brel(q)
    q0_b = (np.log(1/vfactmin-1)+q0*22)/22
    b = b_c/(1+np.exp(-22*(q-q0_b)))
    
    a = a_c
    a = (lcerel>1)*a*fisomrel + (lcerel<=1)*a
    
    #%%
    dvdf_isom_con = b/(q*(fisomrel*fmax+a)) # slope in the isometric point at wrt concentric part
    dvdf_isom_ecc = dvdf_isom_con/slopfac # slope in the isometric point at wrt eccentric part
    dFdvcon0      = 1/dvdf_isom_con
    s_as          = 1/sloplin
    p1 = -(fisomrel*q*(fasymp*fmax - fmax))/(s_as - dFdvcon0*slopfac) 
    p2 =  (fisomrel**2*q**2*(fasymp*fmax - fmax)**2)/(s_as - dFdvcon0*slopfac)
    p3 =  -fasymp*fisomrel*q*fmax;
    p4 =  -s_as
    
    #%% Regions
    r_c1 = ((vce<=0) * (dvdf_isom_con<=sloplin)) # Concentric, dvdf_isom_con<=sloplin (normal)	
    r_c2 = ((vce<=0) * (dvdf_isom_con>sloplin)) # Concentric dvdf_isom_con>sloplin (defective case)
    r_e1 = ((vce>0) * (dvdf_isom_ecc<=(sloplin/slopfac))) # Eccentric, dvdf_isom_ecc<=sloplin (normal) 
    r_e2 = ((vce>0) * (dvdf_isom_ecc>(sloplin/slopfac))) # Eccentric, dvdf_isom_ecc>sloplin (defective case)
    
    #%% Compute CE force
    fce_c1 = (q*(b*fisomrel*fmax+a*vce)) / (b-vce) # Concentric, dvdf_isom_con<=sloplin (normal)
    fce_c2 = q*fisomrel*fmax+vce/sloplin # Concentric dvdf_isom_con>sloplin (defective case)
    fce_e1 = (p2-(p3+p4*vce)*(p1+vce))/(p1+vce) # Eccentric, dvdf_isom_con<=sloplin (normal)
    fce_e2 = q*fisomrel*fmax+((vce*slopfac)/sloplin) # Eccentric, dvdf_isom_con>sloplin (defective case)

    #%% Output
    fce = r_c1*fce_c1 + r_c2*fce_c2 + r_e1*fce_e1 + r_e2*fce_e2
    if 'fmax' in muspar:
        fcerel = fce/muspar['fmax']
    else:
        fcerel = None
    
    return fce, fcerel, [fce_c1, fce_c2, fce_e1, fce_e2], [r_c1, r_c2, r_e1, r_e2]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def ForceEQ(lmtc,gamma,muspar):
    """
    ForceEQ Finds relative CE length such that SEE force equals the sum of CE
        and PEE force, for any given MTC-length and gamma (i.e., normalised
        concentration Ca2+ between the filaments)
    
    Inputs:
        lmtc        =   the muscle-tendon-complex length [m]
        gamma       =   the normalised concentration Ca2+ between the filaments [ ]
        muspar      =   dict with muscle parameter values
    
    Outputs:
        fsee        =   the SEE force [N]
        lcerel      =   the relative CE length, i.e. CE length divided by
                            optimum CE length) [ ]
        fce         =   the CE force [N]
        fpee        =   the PEE force [N]
    """
    
    #%% Unravel parameter values
    lce_opt = muspar['lce_opt']     # [m] CE optimum length (i.e., CE length at which isometric CE force is maximal)
    lsee0   = muspar['lsee0']       # [m] SEE slack length
    fmax    = muspar['fmax']        # [N] maximal CE force
    C       = -1/muspar['w']**2     # [ ] parameter based on width of 2nd order polynomial of isometric CE force-length relation
    ksee    = muspar['ksee']        # [N/m^2]   shape parameter that scales the stiffness of PEE
    
    #%% Get initial guess of relative CE length   
    q = (gamma>muspar['q0'])*ActState(gamma,1,muspar)[0] + (gamma<=muspar['q0'])*muspar['q0'] # [ ]
    # lcerel =  -(ksee*lce_opt*lsee0 - ksee*lce_opt*lmtc + fmax**(1/2)*q**(1/2)*(ksee*C*lce_opt**2 + ksee*C*lmtc**2 + ksee*C*lsee0**2 + ksee*lce_opt**2 - C*fmax*q - 2*ksee*C*lce_opt*lmtc + 2*ksee*C*lce_opt*lsee0 - 2*ksee*C*lmtc*lsee0)**(1/2) + C*fmax*q)/(ksee*lce_opt**2 - C*fmax*q)  # [ ]
    a = ksee*lce_opt*lsee0 - ksee*lce_opt*lmtc + C*fmax*q
    b1 = fmax**(1/2)*q**(1/2)
    b2 = ksee*C*lce_opt**2 + ksee*C*lmtc**2 + ksee*C*lsee0**2 + ksee*lce_opt**2 - C*fmax*q - 2*ksee*C*lce_opt*lmtc + 2*ksee*C*lce_opt*lsee0 - 2*ksee*C*lmtc*lsee0
    b2 = (b2>=0)*b2 + (b2<0)*0   
    b = b1*b2**(1/2) 
    c = ksee*lce_opt**2 - C*fmax*q
    lcerel = -(a+b)/c
    
    #%% Newton root finding to find relative CE length
    # Set values for first round
    fce = 101; fpee = 0; fsee = 0; dlcerel = 0
    # Set tolerance
    tolF = 1e-4; iRound = 0;
    
    # Go
    while (np.max(np.abs(fce+fpee-fsee))>tolF):
        lcerel = lcerel+dlcerel  # [ ]
        fisomrel,Kcerel = ForceLength(lcerel,muspar)  # [ , ]
        q,Kq = ActState(gamma,lcerel,muspar)[0:2] # [ , ]
        fce = q*fisomrel*muspar['fmax'] # [N]
        Kce = fisomrel*Kq*muspar['fmax']+q*Kcerel*muspar['fmax'] # [N]
        lce = lcerel*muspar['lce_opt']
        lsee = lmtc-lce
        lpee = lce
        fsee,fpee,Ksee,Kpee = LEE2Force(lsee,lpee,muspar)[0:4] # [N, N, N/m, N/m]
        Ksee = -Ksee*muspar['lce_opt'] # [N]
        Kpee = Kpee*muspar['lce_opt'] # [N]
        dlcerel = (fce+fpee-fsee)/(Ksee-Kce-Kpee) # [ ]
        iRound = iRound+1
        pass
    
    #%% Output
    return fsee,lcerel,fce,fpee

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def ForceLength(lcerel,muspar):
    """
    ForceLength Computes the relative isometric CE force based on the relative
        CE length.
    
    Inputs:
        lcerel      =   the relative CE length (i.e. CE length divided by 
                            optimum CE length [ ]
        muspar      =   dict with muscle parameter values
    
    Outputs:
        fisomrel    =   the relative isometric CE force (i.e., CE isometric 
                            force normalised by the maximal CE isometric 
                            force) [ ]
        kce         =   dfisomrel/dlcerel [ ]
    """
        
    #%% Unravel parameter values
    n = muspar['n'] # [ ]
    C = -1/muspar['w']**n # [ ]
    
    #%% Compute tails of parabola (exponential tails)
    Fp = 0.1 # exp function kicks in a 10% fisomrel
    xp = ((Fp-1)/C)**(1/2) # intercept of function with y=Fp
    xp = np.array([-xp+1, xp+1]) # two solutions because of root
    dFdx = 2*C*(xp-1) # first derivative of F at xp
    # function has the form: y=a*exp(b*x), so:
    b = dFdx/Fp
    a = Fp/(np.exp(b*xp))

    #%% Compute fisom_rel
    fisomrel = C*(lcerel-1)**n+1 # [ ]
    fisomrel_t1 = a[0]*np.exp(b[0]*lcerel) # [ ] exponential tail for lcerel < xp[0]
    fisomrel_t2 = a[1]*np.exp(b[1]*lcerel) # [ ] exponential tail for lcerel > xp[1]
    fisomrel = (lcerel>=xp[0])*(lcerel<=xp[1])*fisomrel + (lcerel<xp[0])*fisomrel_t1 + (lcerel>xp[1])*fisomrel_t2
    
    #%% Compute derivative (i.e., fisomrel/dlcerel)
    kce = 2*C*(lcerel-1) # [ ] dfisomrel/dlcerel
    kce_t1 = a[0]*b[0]*np.exp(b[0]*lcerel) # [ ] dfisomrel/dlcerel for lcerel < xp[0]
    kce_t2 = a[1]*b[1]*np.exp(b[1]*lcerel) # [ ] dfisomrel/dlcerel for lcerel < xp[0]
    kce = (lcerel>=xp[0])*(lcerel<=xp[1])*kce + (lcerel<xp[0])*kce_t1 + (lcerel>xp[1])*kce_t2  
    
    #%% Output
    return fisomrel,kce

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def GetForceCE(lmtc,fsee,muspar):
    """
    GetForceCE Computes CE force based on MTC length and SEE force.
    
    Inputs:
        fsee        =   SEE force [N]
        fpee        =   PEE force [N]
                
    Outputs:
        lsee        =   SEE length [m]
        lpee        =   PEE length [m]
    """
    
    #%% Unravel parameter values
    lpee0       = muspar['lpee0']       # [m]       PEE slack length
    lsee0       = muspar['lsee0']       # [m]       SEE slack length 
    ksee        = muspar['ksee']        # [N/m^2]   shape parameter that scales the stiffness of PEE
    kpee        = muspar['kpee']        # [N/m^2]   shape parameter that scales the stiffness of SEE
    
    #%% Computations
    # SEE
    esee    = (fsee/ksee)**(1/2)                    # [m] SEE elongation
    lsee    = esee+lsee0                            # [m] SEE length
    
    # 
    lce     = lmtc-lsee                             # [m] CE length
    lpee    = lce                                   # [m] PEE length
    epee    = lpee-lpee0                            # [m] PEE elongation
    fpee    = (epee<0)*0 + (epee>0)*(kpee*epee**2)  # [N] PEE force
    fce     = fsee-fpee                             # [N] CE force
                                  
    #%% Output
    return fce,fsee,fpee,lce,lsee,lpee

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def LEE2Force(lsee,lpee,muspar):
    """
    LEE2Force Computes SEE and PEE force based on SEE and PEE length.
    
    Inputs:
        lsee        =   SEE length [m]
        lpee        =   PEE length [m]
        muspar      =   dict with muscle parameter values
        
    Outputs:
        fsee        =   SEE force [N]
        fpee        =   PEE force [N]
        Ksee        =   dFsee/dLsee [N/m]
        Kpee        =   dFpee/dLpee [N/m]
    """
    
    #%% Unravel parameter values
    lpee0       = muspar['lpee0']       # [m]       PEE slack length
    lsee0       = muspar['lsee0']       # [m]       SEE slack length    
    ksee        = muspar['ksee']        # [N/m^2]   shape parameter that scales the stiffness of PEE
    kpee        = muspar['kpee']        # [N/m^2]   shape parameter that scales the stiffness of SEE
    
    #%% Computations
    # SEE
    esee = lsee-lsee0                               # [m] SEE elongation
    fsee = (esee<0)*0 + (esee>=0)*(ksee*esee**2)    # [N] SEE force
    Ksee = (esee<0)*0 + (esee>=0)*(2*ksee*esee)     # [N/m] dFsee/dLsee
    
    # PEE
    epee = lpee-lpee0                               # [m] PEE elongation
    fpee = (epee<0)*0 + (epee>=0)*(kpee*epee**2)    # [N] PEE force
    Kpee = (epee<0)*0 + (epee>=0)*(2*kpee*epee)     # [N/m] dFpee/dLpee
    
    #%% Output
    return fsee,fpee,Ksee,Kpee