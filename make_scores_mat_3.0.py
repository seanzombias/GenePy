# -*- coding: utf-8 -*-
"""
Created on Wed Sept 27 12:02:04 2017

@author: Enrico Mossotto
"""
import sys, requests
import numpy as np
import collections as cl
from scipy.stats import mannwhitneyu as rks



server = "http://rest.ensembl.org/variation/human/" #Ensembl API server

#%% INPUT files

def import_data():
    IBD=([line.rstrip().split('\t') for line in open(sys.argv[1],'r')])
    LEV=([line.rstrip().split('\t') for line in open(sys.argv[2],'r')])
    return IBD, LEV

cases, levin = import_data()
gene=sys.argv[3]


cases_header=np.array(cases.pop(0))
levin_header=np.array(levin.pop(0))

#%% consolidate databases
### IBD CD CASES  ###
cases=np.array(cases)

scores=cases[:,6:22]
scores_names=cases_header[6:22]

rsids = cases[:,24]
#%%
known_fa = cases[:,22].astype('float')

for i,_x in enumerate(known_fa):
    if np.isnan(_x):
        if np.isnan(cases[i,23].astype(float)):
            continue
        else:
            known_fa[i]=cases[i,23].astype(float)
#%%
proxies = {'http' : "socks5://localhost:6789"} #on if behind a firewall
#proxies = {}

# Obtain frequencies not collected by Annovar using the Ensembl API
freqs = np.zeros((rsids.shape[0], 2))
freqs[:] = np.nan
freqs[:,1]= known_fa
freqs[:,0]= 1-known_fa
def get_alts_freq(r, a, rs):
    if rs != 'nan': # is variant associated to rsid?
        ref = np.nan
        alt = np.nan
        ext = rs+"?pops=1"
#        v = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        v = requests.get(server+ext, headers={ "Content-Type" : "application/json"}, proxies=proxies) #if behind a firewall
        v = v.json()
        fr = 0
        fa = 0
        if '1000Genomes' in v['evidence']:
            for x in v['populations']:
                if x['population']=='1000GENOMES:EUR' and x['allele']==r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population']=='1000GENOMES:EUR' and x['allele']==a and fa == 0:
                    alt = x['frequency']
                    fa = 1
                
        elif '1000Genomes' in v['evidence'] and (fa == 0 or fr == 0):
            for x in v['populations']:
                if x['population'] == '1000GENOMES:phase_3:EUR' and x['allele'] == r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population'] == '1000GENOMES:phase_3:EUR' and x['allele'] == a and fa == 0:
                    alt = x['frequency']
                    fa = 1
                    
        elif 'HapMap' in v['evidence'] and (fa == 0 or fr == 0):
            for x in v['populations']:
                if x['population'] == 'CSHL-HAPMAP:HapMap-CEU' and x['allele'] == r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population'] == 'CSHL-HAPMAP:HapMap-CEU' and x['allele'] == a and fa == 0:
                    alt = x['frequency']
                    fa = 1         
                    
        elif 'ESP' in v['evidence'] and (fa == 0 or fr == 0):
            for x in v['populations']:
                if x['population'] == 'ESP6500:European_American' and x['allele'] == r and fr == 0:
                    ref = x['frequency']
                    fr = 1
                if x['population'] == 'ESP6500:European_American' and x['allele'] == a and fa == 0:
                    alt = x['frequency']
                    fa = 1
    
        if fa == 0 and fr == 0: # add novel freq if no annotation
            ref = 0.99999
            alt = 0.00001
            
        elif fa == 0 and fr == 1:
            alt = 1-fr
            
        elif fa == 1 and fr == 0:
            ref = 1-fa
    else:# add novel freq if no rsid entry
        ref =0.99999
        alt = 0.00001

    return ref, alt

# run the frequency function for each rsid
for i, n in enumerate(rsids):
    _ref, _alt = get_alts_freq(cases[i,2], cases[i,3], n)
    if np.isnan(freqs[i,0]):
        freqs[i,0] = _ref
    if np.isnan(freqs[i,1]):
        freqs[i,1] = _alt
    
freqs[freqs==0]=0.00001 # set 0 MAFs as minimum otherwise log fails


cases=cases[:,25:].astype(float)/2.00 #recode genotypes as 0,0.5,1
cases_header=cases_header[25:]        
#
#
#####  LEVIN CONTROLS ###
levin=np.array(levin)
levin=levin[:,25:].astype(float)/2.00
levin_header=levin_header[25:]

final_sample_list=np.hstack((cases_header,levin_header))

#%% Clalculate the GenePy scores and test

def score_db(ibd,lev,score,freq,sname):
    S=np.copy(score)
    db1=np.copy(ibd)
    db2=np.copy(lev)

    out1=[]
    for i in range(db1.shape[0]):
        if S[i]!='nan':
            deleter=float(S[i])
            db1[i][db1[i]==0.5]=deleter*-np.log10(float(freq[i,0])*float(freq[i,1]))
            db1[i][db1[i]==1]=deleter*-np.log10(float(freq[i,1])*float(freq[i,1]))
            out1.append(db1[i])
    out2=[]
    for i in range(db2.shape[0]):
        if S[i]!='nan':
            deleter=float(S[i])
            db2[i][db2[i]==0.5]=deleter*-np.log10(float(freq[i,0])*float(freq[i,1]))
            db2[i][db2[i]==1]=deleter*-np.log10(float(freq[i,1])*float(freq[i,1]))
            out2.append(db2[i])
    
    out1=np.array(out1)
    out2=np.array(out2)

    
    out1 = np.sum(out1,axis=0)
    out2 = np.sum(out2,axis=0)
    
    tmp = np.hstack((out1,out2))
    gg  =np.array([gene]*len(final_sample_list))
    U = np.vstack((final_sample_list,tmp,gg)).T
    
    return rks(out1,out2), U


print '#Gene\tScore\tIBD_vs_Levin'

for i in range(scores.shape[1]):
    if cl.Counter(scores[:,i])['nan'] < (scores.shape[0]-1): #compute metascores if at least 1 variant
        pvs, U =score_db(cases,levin,scores[:,i],freqs, scores_names[i])
        print '%s\t%s\t%.5f' %(gene,scores_names[i],pvs[1])
        np.savetxt('./'+scores_names[i]+'/'+gene+'_'+scores_names[i]+'_matrix',U, fmt='%s', delimiter='\t')
    else:
        print '%s\t%s\t%s\n' %(gene,scores_names[i],'-')
