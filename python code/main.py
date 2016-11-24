# -*- coding: utf-8 -*-
"""
Created on Tue Oct 08 09:52:03 2013

@author: sophie kolbe
"""

import plot_specs as pl
import identify as id
import score as s
import numpy as np



# +++++++++  MAIN METHODE  +++++++++++

#******************* files   ******************
#file1 = 'random_spec.mgf'
#file2 = 'x_link_charge_2.mgf'
#file3 = 'glufib.mgf'

file1 = 'mtaq_scan_1756.mgf' #MGLPPLLSLPSNAAPR 
file2 = 'mtaq_scan_1164.mgf' #TSVYYLGEVFPQK
file3 = 'mtaq_scan_1036.mgf' #VLEPACAHGPFLR
file4 = 'mtaq_scan_980.mgf' #VLEPACAHGPFLR

#***************  parameters ************
#peptide='EGVNDNEEGFFSAR' #glufib
#peptide='TSVYYLGEVFPQK' #mtaq
peptide='MGLPPLLSLPSNAAPR' #mtaq
#peptide='VLEPACAHGPFLR' #mtaq
peptide_sequence=['VETPPEVVDFMVSLAEAPR','LPDSSLVQWLNSEAMQK','LEISGMPLGDLFHIR',
'MGLPPLLSLPSNAAPR','DFYATPHLVVAHTK','TSVYYLGEVFPQK','VLEPACAHGPFLR',
'EYGFHTSPESAR','NLKPGWVDYEK','EPGPGLVPVLTGR']
#charge
max_charge=1

# modifications
cc1 =  669.23293 #fragmentation 1 of cc: 1116-447
cc2 =  832.35377 #fragmentation 2 of cc: 1115.48922-284.14272
cc3 = 1115.48922 # native cc without fragmentation: 1116 geblitzt
#cc4 = 1143.49537 # native cc without fragmentation: 1143 ungeblitzt - 
#cc5 = 866.42213

mass_list=[0,cc1,cc2,cc3]
modi_list=['cc']#'-amino','-water

#*********************** examples **************************
##
#for file in [file2]:
#    for c in np.arange(1,5): 
#        s.peptide_score(file, [peptide],c,mass_list, '-o')
#         #pl.identify_plot(file,peptide, c,mass_list,1) #theo specs plots of peptide
#for file in [file2]:
#    for c in np.arange(1,5): 
#        print 'charge:',c
#        for mass in mass_list:
#            print 'mass:',mass
#            s.xcorr(file,peptide,c,['cc'],mass,'-o') 

 
# peptide list       
#for peptide in ['VLEPACAHGPFLR']:  
#    print ''
#    print 'peptide => ',peptide   
#    #file list
#    for file in [file3]:   
#        print '     file::',file
#        #charge 1-4
#        for c in np.arange(1,5):
#            print '             charge=',c
#            pl.identify_plot(file,peptide, c,mass_list,1)
#            s.sequence_score(file, peptide,c,modi_list,mass_list,'-p')
#            for mass in [0,cc1,cc2,cc3]:
#            #mass list
#                print '                 mass ',mass
#                id.identify(file,peptide,c,['cc'],mass,'-o')
#                print ' '
#            print ''
#        print ''
#        print ''



