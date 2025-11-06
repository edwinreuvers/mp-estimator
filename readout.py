# -*- coding: utf-8 -*-

#%% Load packages
import os, pickle
import pandas as pd

# Set directories
cwd = os.path.dirname(os.path.abspath(__file__))
baseDir = os.path.join(cwd,'')
dataDir = os.path.join(baseDir,'data')
parDir = os.path.join(baseDir,'parameters')

#%% Example 1: Estimated CE, SEE and PEE force-length parameter values only
results = {}
# Loop through and load each dict
for mus in ['GMs1', 'GMs2', 'GMs3']:
    filepath = os.path.join(parDir,mus,f'{mus}_example1.pkl')
    
    estpar = pickle.load(open(filepath, 'rb'))
    results[mus] = estpar
df_ex1 = pd.DataFrame.from_dict(results, orient='index')
print(df_ex1)

#%% Example 2: Estimated CE force-velocity parameter values only
results = {}
# Loop through and load each dict
for mus in ['GMs1', 'GMs2', 'GMs3']:
    filepath = os.path.join(parDir,mus,f'{mus}_example2.pkl')
    
    estpar = pickle.load(open(filepath, 'rb'))
    results[mus] = estpar
df_ex2 = pd.DataFrame.from_dict(results, orient='index')
print(df_ex2)

#%% Example 3: Estimated contraction dynamics parameter values using the 'traditional method'
results = {}
# Loop through and load each dict
for mus in ['GMs1', 'GMs2', 'GMs3']:
    filepath = os.path.join(parDir,mus,f'{mus}_example3.pkl')
    
    estpar = pickle.load(open(filepath, 'rb'))
    results[mus] = estpar
df_ex3 = pd.DataFrame.from_dict(results, orient='index')
print(df_ex3)

#%% Example 4: Estimated contraction dynamics parameter values using the 'improved method'
results = {}
# Loop through and load each dict
for mus in ['GMs1', 'GMs2', 'GMs3']:
    filepath = os.path.join(parDir,mus,f'{mus}_example4.pkl')
    
    estpar = pickle.load(open(filepath, 'rb'))
    results[mus] = estpar
df_ex4 = pd.DataFrame.from_dict(results, orient='index')
print(df_ex4)

#%% Example 5: Estimated excitation dynamics parameter values 
# with contraction dynamics parameter using the the 'improved method'
results = {}
# Loop through and load each dict
for mus in ['GMs1', 'GMs2', 'GMs3']:
    filepath = os.path.join(parDir,mus,f'{mus}_example5_IM.pkl')
    
    estpar = pickle.load(open(filepath, 'rb'))
    results[mus] = estpar
df_ex5 = pd.DataFrame.from_dict(results, orient='index')
print(df_ex5)