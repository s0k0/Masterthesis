# -*- coding: utf-8 -*-
"""
Created on Wed Nov 06 10:43:04 2013

@author: sophie kolbe
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 12:44:52 2013

@author: sophie kolbe
"""
from pyteomics import mass
from pyteomics import mgf

import numpy as np

#read experimental spectrum from mgf-format file; only single spectras accepted!
def read_mgf(file):  
    spectra = mgf.read(file)      
    spectrum = next(spectra)
    #create parameters: mz = m/z array, para = parameters, amp = intensity array
    mz = spectrum['m/z array']
    params= spectrum['params']
    amp = spectrum['intensity array']
    return(mz,params,amp) 


#creates theoretical spectras of y-ions series fom given peptide  
def fragments_y(peptide, maxcharge=1):
    for i in xrange(0, len(peptide)): #changed  to catch ending aas
        yield mass.fast_mass(peptide[i:], ion_type='y', charge=maxcharge)
        
#creates theoretical spectras of b-ions series fom given peptide                            
def fragments_b(peptide, maxcharge=1):
    for i in xrange(1, len(peptide)+1): #changed  to catch ending aas
        yield mass.fast_mass(peptide[:i], ion_type='b', charge=maxcharge)
        

#compares theoretical spectra to experimental spectra and detect hits, error threshold defines the sensitivity.
def compare2theo(theo_spectra, mz_values,amp_values,error_threshold):
    match_mz=[]
    match_amp=[]
    match_hits=[]
    #scan theo spectra for matches in experimental spectra
    for theo_peak in theo_spectra:
        #get nearest match in mz values
        closest= min( mz_values, key=lambda x:abs(x-theo_peak))
        #if error threshold is not crossed, append hit, mz and amp to list
        if abs(closest-theo_peak)<=error_threshold:
            match_hits.append(theo_peak)
            match_mz.append(closest)
            match_amp.append(amp_values[mz_values.tolist().index(closest)])

    return (match_hits,match_mz,match_amp)

#reports matches between given peptide and experimental spectra form file. error threshold can be varried.
#depending on given modi (list [] containing 'cc','-water' and/or '-amino'), different ion series are tracked.
# with the command '-o' the results are directly reportet, else returned in variables y_seq,y_amp or y_mz and
#b_seq,b_amp,b_mz
def identify(file, peptide,charge, modi_type,mass, command='', error=0.5):   
    #read data from given file    :
    mz,params,amplitudes= read_mgf(file)
    # normalize amplitudes and display them as relativ amplitudes in %
    amp = (amplitudes/max(amplitudes))*100
    
    if modi_type==['cc'] and mass==0:
        # 1. create theoretical spectras of peptide sequence to given charge
        theo_y= list(fragments_y(peptide,maxcharge=charge))
        theo_b= list(fragments_b(peptide,maxcharge=charge))
   
        # 2. compare theoretical to experimental spectrum and return matching peaks 
        hit_y,match_mz_y,match_amp_y = compare2theo(theo_y,mz,amp,error)
        hit_b,match_mz_b,match_amp_b = compare2theo(theo_b,mz,amp,error)
        label = 'unmodified'        
    
    else:   
        diff_mass=0
        label=''
        if 'cc' in modi_type:
            #add mass of capture compound (cc) ions
            diff_mass += mass/charge
            if mass==0:
                 label = label + 'native'  
            else:
                label = label+ 'cc fragment:' + str(int(mass))

        if '-water' in modi_type:
            #subtract mass of neutral loss of water ions
            diff_mass += -18.01057/charge
            label = label + ' -water '

        if '-amino' in modi_type:
            #subtract mass of neutral loss of amino ions
            diff_mass += -17.02655/charge
            label = label + ' -amino '

         
        t_y= list(fragments_y(peptide,maxcharge=charge))
        t_b= list(fragments_b(peptide,maxcharge=charge))
        
        #alter theoretical spectra according to mass modifications
        theo_y = [x+diff_mass for x in t_y]
        theo_b = [x+diff_mass for x in t_b]   

        hit_y,match_mz_y,match_amp_y = compare2theo(theo_y,mz,amp,error)
        hit_b,match_mz_b,match_amp_b = compare2theo(theo_b,mz,amp,error)
 
    y_seq=[]  
    y_mz=[]
    y_amp=[]
  
    #fill results of compare2theo into peptide scheme and add gaps when neccessairy
    for position in np.arange(0,len(peptide)):
            if theo_y[position] in hit_y:
                index = hit_y.index(theo_y[position])
                y_seq.append(peptide[position])
                y_mz.append(match_mz_y[index])
                y_amp.append(match_amp_y[index])  
            else:
                y_seq.append('-')
                y_mz.append(0)
                y_amp.append(0)

    
    b_seq=[]
    b_mz=[]  
    b_amp=[]
    #fill results of compare2theo into peptide scheme and add gaps when neccessairy
    for position in np.arange(0,len(peptide)):
            if theo_b[position] in hit_b:
                index = hit_b.index(theo_b[position])
                b_seq.append(peptide[position])
                b_mz.append(match_mz_b[index])
                b_amp.append(match_amp_b[index])  
            else:
                b_seq.append('-')
                b_amp.append(0)
                b_mz.append(0)

    #output
    if command=='-o':
        print label
        print 'y:',y_seq
        print 'amp:',y_amp
        print 'm/z:',y_mz
        print 'b:',b_seq
        print 'amp:',b_amp
        print 'm/z:',b_mz
    return(y_seq,y_mz,y_amp,b_seq,b_mz,b_amp)
