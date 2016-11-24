# -*- coding: utf-8 -*-
"""
Created on Tue Oct 01 11:38:49 2013

@author: sophie kolbe
"""


#import matplotlib as plt
import pylab as py
import numpy as np
import identify as id
import matplotlib.pyplot as plt
   

#plots the results of identify()  
def identify_plot(file, peptide,charge, modi_mass=[], zoom=1):

     #read file 
     values, params, amplitude= id.read_mgf(file)
     #procentual amplitdues for relative abundance
     ampl = (amplitude/max(amplitude))*100 
     #offset and scales for zoom into the spectra
     if zoom>0 and zoom<=1:
         offset = (max(ampl)*zoom)/4
         y_lim_scale = max(ampl)*zoom
         one_per = y_lim_scale/100
     # create theoretical fragmentation of given peptide
     theo_y= list(id.fragments_y(peptide,maxcharge=charge))
     theo_b= list(id.fragments_b(peptide,maxcharge=charge))
         
     # create plot title
     headline = 'File: '+ file+ ', Peptide: '+ peptide + ' (charge '+ str(charge)+')'
     
     #headline = 'File: '+ file
     #create figure
     py.figure()
     py.xlabel('m/z')
     py.ylabel('Intensity, rel. units %')
     py.title(headline)
     
     # plot experimental spectra in grey
     py.bar(values, ampl, width=2.5, linewidth=0.1, color='grey',alpha=0.3, label='data')
     
     # some label and axis configs
     py.xlim(xmin=0)
     py.ylim(ymax= y_lim_scale)
     py.xticks(np.arange(0, max(values)+1, 50))
     py.tick_params(axis='both', which='major', labelsize=8)
     py.tick_params(axis='both', which='minor', labelsize=5)
#     
     
     # --- 0. benchmarks
     #0.1 plot theoretical spectras if no modi given
     if modi_mass==[]:
         py.bar(theo_y,np.ones(len(theo_y))*offset, linewidth=0.1, width=2.5, alpha=0.3, color='red',label='theo y')
         py.bar(theo_b,np.ones(len(theo_y))*offset, linewidth=0.1, width=3.5, alpha=0.3 ,color='lightgreen', label='theo b')
    
        
        # 0.2 annotate peaks
         #0.2.1 theo y ions
         for i in np.arange(0,len(peptide)):     
             py.annotate(peptide[i], xy=(theo_y[i], offset), 
                         xytext=(theo_y[i], offset +2*one_per),
                         bbox=dict(boxstyle='round,pad=0.2', fc='red', alpha=0.5 ))
         # 0.2.1 theo b ion            
             py.annotate(peptide[i], xy=(theo_b[i],offset), 
                         xytext=(theo_b[i], offset +6*one_per),
                         bbox=dict(boxstyle='round,pad=0.2', fc='lightgreen', alpha=0.5 ))
         
      
     alpha = 0.75
     for mass in modi_mass:              
         # --- 1. mono matches 
         if mass==0:   
             # 1.0 use def identify to find matching aa sequences 
             y_seq,y_mz,y_amp,b_seq,b_mz,b_amp = id.identify(file,peptide,charge, ['cc'], 0) 
             
            # 1.1 plot peaks
             py.bar(y_mz,y_amp, width=5, linewidth=0.1, color='red', label='y [M+1H]+'+str(charge))
             py.bar(b_mz,b_amp, width=5,linewidth=0.1 ,color='lightgreen', label='b [M+1H]+'+str(charge))
            
             # 1.2  annotate labels to peak hits
             # 1.2.1 y ion label
    
             for i in np.arange(0,len(peptide)):
                 if y_seq[i]=='-':
                     py.annotate(y_seq[i], xy=(theo_y[i], 0), 
                             xytext=(theo_y[i], offset +4*one_per))
                 else:
                     py.annotate(y_seq[i], xy=(y_mz[i], offset), 
                             xytext=(y_mz[i], offset +4*one_per),
                             bbox=dict(boxstyle='round,pad=0.2', fc='red' ))
        
                             
             # 1.2.2 b ion label 
                 if b_seq[i]=='-':
                     py.annotate(b_seq[i], xy=(theo_b[i], 0), 
                             xytext=(theo_b[i], offset +8*one_per))
                 else:
                     py.annotate(b_seq[i], xy=(b_mz[i], offset), 
                             xytext=(b_mz[i], offset +8*one_per),
                             bbox=dict(boxstyle='round,pad=0.2', fc='lightgreen'))
        
     # --- 2. cc matches   
         else:   
             #compute sequence
             y_seq, y_mz,y_amp, b_seq,b_mz,b_amp =id.identify(file,peptide,charge,['cc'],mass)
             # 2.1. plot peaks
             py.bar(y_mz,y_amp, width=5, linewidth=0.1, color='orange', label=  'y + '+str(int(mass))+'[M+1H]+'+str(charge), alpha=alpha)
             py.bar(b_mz,b_amp, width=5, linewidth=0.1 ,color='blue', label='b + '+str(int(mass) )+'[M+1H]+'+str(charge),alpha=alpha)
            
            # 2.2.annotate labels to peak 
             # 2.2.1 y ions
             offset = offset + 12*one_per
             for i in np.arange(0,len(peptide)):
                 if y_seq[i]=='-':
                     py.annotate(y_seq[i], xy=(theo_y[i]+(mass/charge), 0), 
                         xytext=(theo_y[i]+(mass/charge), offset))
                 else:
                     py.annotate(y_seq[i], xy=(y_mz[i],offset), 
                         xytext=(y_mz[i],offset),
                         bbox=dict(boxstyle='round,pad=0.2', fc='orange', alpha=alpha))
            
            # 2.2.2 b ions            
                 if b_seq[i]=='-':
                     py.annotate(b_seq[i], xy=(theo_b[i]+(mass/charge), 0), 
                         xytext=(theo_b[i]+(mass/charge), offset +4*one_per))
                 else:
                     py.annotate(b_seq[i], xy=(b_mz[i],offset), 
                             xytext=(b_mz[i], offset+4*one_per ),
                             bbox=dict(boxstyle='round,pad=0.2', fc='blue', alpha=alpha))            
             alpha = alpha - 0.25  

     #--- 3. reporter matches
     reporter=[284,447]
     r_val=[]
     r_amp=[]
     # 3.1 check for occurence
     for i in reporter:
         closest= min(values, key=lambda x:abs(x-i))
         if round(closest)== i:
             r_val.append(closest)
             r_amp.append(ampl[values.tolist().index(closest)])
     # 3.2 plot occurence
     if r_val != []:
         py.bar(r_val,r_amp,width=5,color='yellow', label='reporter')
     # 3.2.1 annotate occurence
         for i in np.arange(0,len(r_val)):
             py.annotate(int(round(r_val[i])), xy=(r_val[i],r_amp[i]), ha='center',
                         xytext=(r_val[i],r_amp[i]+2*one_per),
                         bbox=dict(boxstyle='round,pad=0.2', fc='yellow')) 
    
     py.legend(loc=2 , prop={'size':11})
     py.show()       
         
   
#plot results of sequence_score()
def sequence_plot(file,peptide,charge,unmodified, modified,compare):
    
     #write index for peptides on x axis     
     x_index=np.arange(0,len(peptide))
     #create figure
     plt.figure()
    
     plt.suptitle(file+', '+peptide+', until charge:'+str(charge))
     #set ticks to peptide sequence
    
     plt.subplot(211)
     plt.title('Amplitudes of Hits')
     plt.xlabel('Position') 
     plt.xticks(x_index,peptide)
     plt.ylabel('rel. Occurence in %')
     plt.plot(x_index,unmodified[0], 'r',  linewidth=2, label='y unmodified')
     plt.plot(x_index,unmodified[1], 'g',linewidth=2, label='b unmodified')   
     plt.plot(x_index,modified[0],  'r--' , linewidth=2,label='y modified')
     plt.plot(x_index,modified[1], 'g--', linewidth=2,label='b modified')
     plt.legend()
     
     plt.subplot(212)
     plt.title('Sum of Differences (modified-unmodifed) in Amplitudes of Hits')
     plt.xlabel('Position') 
     plt.xticks(x_index,peptide)
     plt.ylabel('Cummulated Differences')
     plt.plot(x_index,compare[0], color='r',linewidth=3,  label='diff y')
     plt.plot(x_index,compare[1], color='g',linewidth=3, label='diff b')
     plt.axhline(0)   
     plt.legend(bbox_to_anchor=[1, 0.3])
     plt.show()        

      
          