# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:25:14 2014

@author: sophie kolbe
"""

import pylab as py
import numpy as np
import identify as id
import scipy.sparse as sp
import plot_specs as pl



#shapes results of identify() into a matrix, containing the intensities of matching peaks 
# in one row per ion (b,y). Non-matches receive zero entries.    
def score_distribution(file, peptide, charge, modi_list=['cc'], mass_list=[0]):
    
    peptide_scores=[]
   # print 'peptide: ',peptide,', charge:',charge, ' , mass:', mass_list, 'modi:',modi_list

    for mass in mass_list:
            y_scores=[]
            b_scores=[]
            # fetch matching hits
            y_seq,y_mz,y_amp, b_seq,b_mz,b_amp = id.identify(file,peptide,charge,modi_list,mass)
    
            # loop over given peptide match according to localise()
            for i in np.arange(0,len(peptide)):
                #if not a gap, rise temp score and scale with amplitude of hit
                if y_seq[i]!= '-':
                    temp_score=y_amp[i]
                else:
                    temp_score= 0
                #assemble list of temp scores, building a seprate row in matrix
                y_scores.append(temp_score)
           
        
                if b_seq[i]!= '-':
                    temp_score=b_amp[i]
                else:
                    temp_score = 0
                b_scores.append(temp_score)
            
            #append temp scores to peptide score list
            peptide_scores.append(b_scores)
            peptide_scores.append(y_scores)
    # cast peptide score list of arrays together as matrix
    A = np.matrix(peptide_scores)
    return(A)
    
#reports distribution of matching peaks intensities. Computes matrix of amplitudes for
# zero-modification (native) and cc-modifications (average over all given modi masses).
# the command '-o' displays a plot showing the intensities distribution in all positions of
#  sequence and their progress over the charge states.
def sequence_score(file, peptide, max_charge, modi_list=['cc'], mass_list=[0], command=''):
    print 'File:', file, ', Peptide:',peptide,', maxcharge:', max_charge    
    #dummies with zero entries to make summing up easier
    matrix_cc=np.zeros((2, len(peptide)))
    matrix_0=np.zeros((2, len(peptide)))
    #for every charge 
    for charge in np.arange(1,max_charge+1):
        #and ever mass loop over peptide and record hits with according amplitude (%) in matrix
                for mass in mass_list:
                    if mass==0:
                        #score matrix for unmodified ions
                        matrix_0 += score_distribution(file, peptide, charge, modi_list,[0])
                    else:
                        #score matrix for modified ions
                        # remark: this sums _ALL_ hits of modications in mass list!
                        matrix_cc +=  score_distribution(file,peptide,charge,modi_list,[mass])         
    # if more than one modifcation and non-empty scores, normalize and compute average over all modifications. 
    #this way the results are more comparable to unmodified counts! 

    if sum(sum(matrix_cc))>0 and len(mass_list)>1:
        score_matrix_cc = matrix_cc/(len(mass_list)-1)
 
    compare_matrix=[]
    for i in np.arange(0,2):
        row =[]
        diff_score=0
        for j in np.arange(0,len(peptide)):
            diff_score += round(matrix_cc[i][j] -matrix_0 [i][j],1)
            row.append(diff_score) 
        compare_matrix.append(row)
        
        
    if command =='-o':
        
        print 'unmodified matrix:',matrix_0
        print ''
        print 'modified matrix:',score_matrix_cc    
        print ''
        print 'compare matrix:',compare_matrix
    
    if command=='-p':
        pl.sequence_plot(file,peptide,charge,matrix_0, matrix_cc,compare_matrix)     
    return(compare_matrix,matrix_0, score_matrix_cc)


# computes preliminary score for weighting quality of matches 
#between exprimental data and a list of peptides. command '-o' prints
#every single peptide score of list, '-b' reports only best score.
# anyway both of them returned by the function in variables 'best_score'
# and 'return_list'.
def peptide_score(file, peptide_list,max_charge,mass_list, command=''):
    
    #item holder for output:
    score_list=[]
    # loop over every peptide, charge status, modification type and mass patter
    for peptide in peptide_list:
        peptide_score=0
        for mass in mass_list:
            
            # fetch matching hits
            y_seq,y_mz,y_amp, b_seq,b_mz,b_amp = id.identify(file,peptide,max_charge,['cc'],mass)
    
            # counters
            score=0
            consecutive_hits=0
            number_of_hits=0
           # raise the score through summing the matches for y series, scaled over their amplitude
            for i in np.arange(0,len(peptide)):
                if y_seq[i]!= '-':
                    number_of_hits+=1
                    score +=y_amp[i]
                    # if adjacent series, raise counter
                    if i >0 and y_seq[i-1]!= '-':
                        consecutive_hits +=1
                if b_seq[i]!= '-':
                    number_of_hits+=1
                    score +=b_amp[i] 
                    if i >0 and b_seq[i-1]!= '-':
                        consecutive_hits +=1
                        
            # compute final score: sum up intensities of hits (peptide_score), multiply with number of hits,
        # add bonus for consecutive hits and divide through number of possible hits (len(peptide)= #ions)            
            peptide_score += (score*number_of_hits*(1+consecutive_hits))/len(peptide)
            #print peptide+' at '+str(int(mass))+': '+str(int(peptide_score))+ ' and consecutive hits:'+str(consecutive_hits)    
        
        score_list.append((peptide_score))  
        
        #output single computation results 
        if command=='-o':     
            print 'peptide:',peptide, ' -> score:', int(peptide_score)
     #find max value in scores
    if len(peptide_list)>1:   
        max_score = max(score_list)  
        max_score_peptide= peptide_list[score_list.index(max_score)]
        best_score=[max_score_peptide,max_score]
    else:
        best_score=[peptide_list,score_list]
      
    #report best score  
    if command=='-o' or command=='-b':
        print 'best score: ', best_score
        
    return_list=zip(peptide_list,score_list) 
        
    return(best_score,return_list)  
    



#computes pearson correlation coefficient between synthetic peaks and experimental matches 
#focusing their mz-values. command '-o' return results.
def xcorr(file,peptide,max_charge, modi_list=['cc'], mass=0,command='-o',error=0.5):
    #read file
    mz,params,amplitudes= id.read_mgf(file)
    #create theoretical spectra of peptide
    theo_y= list(id.fragments_y(peptide,maxcharge=max_charge))
    theo_b= list(id.fragments_b(peptide,maxcharge=max_charge))
    #fetch data
    y_seq,y_mz,y_amp, b_seq,b_mz,b_amp = id.identify(file,peptide,max_charge,modi_list,mass,'',error)
  
 
    # correlation coefficient; return matrix with XX,XY,YX,XX.  we need xy, so enter only this part of matrix 
    xcorr_y = np.round(np.corrcoef(theo_y,y_mz)[0,1],3)
    
    xcorr_b= np.round(np.corrcoef(theo_b,b_mz)[0,1],3)
    
    if command=='-o':
        print 'seq: ',peptide, '==> xcorr y:', xcorr_y,', xcorr b:', xcorr_b  
    return(xcorr_y,xcorr_b)
 
  