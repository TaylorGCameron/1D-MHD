# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 15:02:00 2018

@author: Taylor
"""

from __future__ import absolute_import
#from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import numpy as np
from numpy import linalg as LA
import matplotlib.dates as mdate
import datetime
import argparse
import struct
#import json # for pretty output
from ai import cdas


#Goes through some text and finds the value associated with some keyword.
def get_keyword(text, keyword):
    ind1 = text.find(keyword)
    ind2 = text[ind1:].find('\n') + ind1
    return text[ind1:ind2][text[ind1:ind2].find('=')+1:].replace(' ','')

#modifies amrvac.par. values is a list of keywords and values to insert. Form is [[keyword, value],[keyword, value],...
def write_keywords(values):
    parfile = open('amrvac.par','r')
    par_text = parfile.read()
    parfile.close()
    
    for i in range(len(values)):
        keyword = values[i][0]
        value = values[i][1]
        ind1 = par_text.find(keyword)
        ind2 = par_text[ind1:].find('\n') + ind1
        par_text = par_text[0:ind1+len(keyword)] + '= '+ value + par_text[ind2:]
    
    parfile = open('amrvac.par', 'w')
    parfile.write(par_text)
    parfile.close()


def pull_ACE(t1, t2):
    
    swe_data = cdas.get_data('sp_phys', 'AC_H0_SWE', mdate.num2date(t1), mdate.num2date(t2), ['Np', 'Vp', 'Tpr', 'V_GSE', 'SC_pos_GSE'])

    mfi_data = cdas.get_data('sp_phys','AC_H3_MFI',mdate.num2date(t1), mdate.num2date(t2),['BGSEc'])
    
    return swe_data, mfi_data

def average_mfi(swe_data, mfi_data):

    Bx_resampled = np.zeros(len(swe_data['EPOCH']))
    By_resampled = np.zeros(len(swe_data['EPOCH']))
    Bz_resampled = np.zeros(len(swe_data['EPOCH']))

    for i in range(len(swe_data['EPOCH'])):
        Bx_resampled[i] = np.nanmean(np.array(mfi_data['BX_GSE'])[np.logical_and(np.array(mfi_data['EPOCH']) > swe_data['EPOCH'][i] - datetime.timedelta(0,32), np.array(mfi_data['EPOCH']) < swe_data['EPOCH'][i] + datetime.timedelta(0,32))])
        By_resampled[i] = np.nanmean(np.array(mfi_data['BY_GSE'])[np.logical_and(np.array(mfi_data['EPOCH']) > swe_data['EPOCH'][i] - datetime.timedelta(0,32), np.array(mfi_data['EPOCH']) < swe_data['EPOCH'][i] + datetime.timedelta(0,32))])
        Bz_resampled[i] = np.nanmean(np.array(mfi_data['BZ_GSE'])[np.logical_and(np.array(mfi_data['EPOCH']) > swe_data['EPOCH'][i] - datetime.timedelta(0,32), np.array(mfi_data['EPOCH']) < swe_data['EPOCH'][i] + datetime.timedelta(0,32))])

    mfi_data['EPOCH'] = swe_data['EPOCH']
    mfi_data['BX_GSE'] = Bx_resampled
    mfi_data['BY_GSE'] = By_resampled
    mfi_data['BZ_GSE'] = Bz_resampled
    return mfi_data

#Function to combine ACE_MFI and ACE_SWE data. (MFI needs to be at the same time resolution as swe)

def combine_ACE(swe_data, mfi_data):

    if len(swe_data['EPOCH']) != len(mfi_data['EPOCH']):
        print('Resample MFI data first!')
        return -1
    
    dtype=np.dtype([('t', '<f8'), ('pos', '<f8', (3,)), ('Bx', '<f8'), ('By', '<f8'), ('Bz', '<f8'), ('vx', '<f8'), ('vy', '<f8'), ('vz', '<f8'), ('n', '<f8'), ('T', '<f8')])
    ACE = np.ndarray(len(swe_data['EPOCH']), dtype = dtype)
    
    ACE['t'] = mdate.date2num(swe_data['EPOCH'])
    ACE['Bx'] = mfi_data['BX_GSE']
    ACE['By'] = mfi_data['BY_GSE']
    ACE['Bz'] = mfi_data['BZ_GSE']
    ACE['vx'] = swe_data['VX_(GSE)']
    ACE['vy'] = swe_data['VY_(GSE)']
    ACE['vz'] = swe_data['VZ_(GSE)']
    ACE['n'] = swe_data['H_DENSITY']
    ACE['T'] = swe_data['H_TEMP_RADIAL']
    ACE['pos'] = np.transpose([swe_data['ACE_X-GSE'], swe_data['ACE_Y-GSE'], swe_data['ACE_Z-GSE']]) 
    
    return ACE

#Cleans up bad values. Flags them as NANS, and then fills them in
def clean_ACE(ACE):
    #Flag values as NANs    
    ACE['vx'][ACE['vx'] < -10**30] = np.nan
    ACE['vy'][ACE['vy'] < -10**30] = np.nan
    ACE['vz'][ACE['vz'] < -10**30] = np.nan
    ACE['Bx'][ACE['Bx'] < -10**30] = np.nan
    ACE['By'][ACE['By'] < -10**30] = np.nan
    ACE['Bz'][ACE['Bz'] < -10**30] = np.nan
    ACE['n'][ACE['n'] < -10**30] = np.nan
    ACE['T'] [ACE['T']  < -10**30] = np.nan
    
    
    #Fill in NANs with the value around them
    ACE['vx'] = fill_nan_t(ACE['t'], ACE['vx'])
    ACE['vy'] = fill_nan_t(ACE['t'], ACE['vy'])
    ACE['vz'] = fill_nan_t(ACE['t'], ACE['vz'])
    ACE['Bx'] = fill_nan_t(ACE['t'], ACE['Bx'])
    ACE['By'] = fill_nan_t(ACE['t'], ACE['By'])
    ACE['Bz'] = fill_nan_t(ACE['t'], ACE['Bz'])
    ACE['n'] = fill_nan_t(ACE['t'], ACE['n'])
    ACE['T']  = fill_nan_t(ACE['t'], ACE['T'] )
    
    
    
    return ACE
        
#Function used in rotating data

def getMVAB0Normal(t1, t2):
    swe_data, mfi_data = pull_ACE(t1, t2)
    mfi_data = average_mfi(swe_data, mfi_data)
    ACE = combine_ACE(swe_data, mfi_data)
    ACE = clean_ACE(ACE)

    Bx = ACE['Bx']
    By = ACE['By']
    Bz = ACE['Bz']

    #plt.plot(By)

    B = np.transpose(np.array([Bx, By, Bz]))
    M = np.zeros([3,3])
    P = np.zeros([3,3])
    
    B_av = np.array([np.nanmean(B[:,0]), np.nanmean(B[:,1]),np.nanmean(B[:,2])])
    e = B_av/ np.sqrt(B_av[0]**2+B_av[1]**2+B_av[2]**2)

    #Create a covariance matrix
    
    for i in range(3):
        for j in range(3):
            M[i,j] = np.nanmean(B[:,i]*B[:,j]) -np.nanmean(B[:,i])*np.nanmean(B[:,j])
            if i == j:
                P[i,j] = 1 - (e[i]*e[j])
            else:
                P[i,j] = -1* (e[i]*e[j])
                    
    #Get eigenvalues and eigenvectors
    M = np.dot(np.dot(P,M),P)
    eigenvalues, eigenvectors = LA.eig(M)
    args = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[args]
    eigenvectors = eigenvectors[:,args]
    
    #The vector corresponding to the middle (absolute value) eigenvalue is the minimum variance direction 
    
    front_normal = eigenvectors[:,1]
    
    #The x component of the vector should point towards the sun (positive)
    if front_normal[0] < 0:
        front_normal = -1*front_normal
    
    # Do a test. For the result to be valid, the second smallest eigenvalue should be x (5 for now) times larger than the smallest
      
    ratio = eigenvalues[2] / eigenvalues[1]   
    n = front_normal  
#    print 'r is ', ratio
    return n, ratio


def rotate_ACE(ACE, angle, axis):
    
    for i in range(0, len(ACE['t'])):
        p_r = rotate(np.array([ACE['pos'][i,0], ACE['pos'][i,1], ACE['pos'][i,2]]), angle, axis)
        B_r = rotate(np.array([ACE['Bx'][i], ACE['By'][i], ACE['Bz'][i]]), angle, axis)
        v_r = rotate(np.array([ACE['vx'][i], ACE['vy'][i], ACE['vz'][i]]), angle, axis)
        
        #if i == 0:
        #    print(p_r[0])
        #    print(ACE['pos'][0,0])
        #    print self.ACE['vy'][i]
        #    print self.ACE['vz'][i]
        #    print v_r

        ACE['pos'][i] = p_r
        
        ACE['Bx'][i] = B_r[0]
        ACE['By'][i] = B_r[1]
        ACE['Bz'][i] = B_r[2]
    
        ACE['vx'][i] = v_r[0]
        ACE['vy'][i] = v_r[1]
        ACE['vz'][i] = v_r[2]
        
        
    return ACE

def get_axis(v1, v2): 
    return np.cross(v1, v2)

def mag(v):
     return np.sqrt(v.dot(v))

def get_angle(v1, v2):
    return np.arccos((np.dot(v1, v2))/(mag(v1)*mag(v2)))



def rotate(v, angle, axis):
    axis = axis/mag(axis)
    
    R = np.array([[(np.cos(angle)+axis[0]**2*(1-np.cos(angle))), (axis[0]*axis[1]*(1-np.cos(angle)) - axis[2]* np.sin(angle)), (axis[0]*axis[2]*(1-np.cos(angle))+axis[1]*np.sin(angle))],
                   [(axis[1]*axis[0]*(1-np.cos(angle)) + axis[2]*np.sin(angle)), (np.cos(angle) + axis[1]**2*(1-np.cos(angle))), (axis[1]*axis[2]*(1-np.cos(angle))-axis[0]*np.sin(angle))],
                   [(axis[2]*axis[0]*(1-np.cos(angle))-axis[1]*np.sin(angle)), (axis[2]*axis[1]*(1-np.cos(angle))+axis[0]*np.sin(angle)), (np.cos(angle)+axis[2]**2*(1-np.cos(angle)))]])
        
    return np.dot(R, v)

#Pulls OMNI data, and puts it in a useful form
def pull_OMNI_BSN(t1,t2):
    data = cdas.get_data('sp_phys', 'OMNI_HRO_1MIN', mdate.num2date(t1), mdate.num2date(t2), ['BSN_x','BSN_y','BSN_z'])
    dtype=np.dtype([('t', '<f8'), ('bsn_x_gse','<f8'),('bsn_y_gse','<f8'),('bsn_z_gse','<f8')])
    BSN = np.ndarray(len(data['EPOCH_TIME']), dtype = dtype)
    
    BSN['bsn_x_gse'] = data['X_(BSN),_GSE']
    BSN['bsn_y_gse'] = data['Y_(BSN),_GSE']
    BSN['bsn_z_gse'] = data['Z_(BSN),_GSE']
    for i in range(len(BSN)):
        BSN['t'][i] = mdate.date2num(data['EPOCH_TIME'][i])
    return BSN

#Functions for loading 1d amrvac data





#Fill in Nans with interpolated data
    
#This is where we interpolate over nans?

def fill_nan(A):
    '''
    interpolate to fill nan values
    '''
    from scipy import interpolate
    inds = np.arange(A.shape[0])
    good = np.where(np.isfinite(A))
    if len(good[0]) == 0:
        return np.nan_to_num(A)
    f = interpolate.interp1d(inds[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(inds))
    return B

def fill_nan_t(t, A):
    '''
    interpolate to fill nan values
    '''
    from scipy import interpolate
    good = np.where(np.isfinite(A))
    if len(good[0]) == 0:
        return np.nan_to_num(A)
    f = interpolate.interp1d(t[good], A[good], bounds_error=False)
    B = np.where(np.isfinite(A), A, f(t))
    return B

#Functions for MVAB0 related stuff for paper
    
def getMVAB0ShiftedACE(ACE, WIND, t1,t2, parameter, interval = 60):
    Re = 6371.
    
    n1 = np.where(ACE['t'] > t1)[0][0]
    n2 = np.where(ACE['t'] < t2)[0][-1]
    
    #Make some variables
    
    normals = np.zeros([len(range(n1,n2)),3])
    ratios = np.zeros([len(range(n1,n2))])
    #print('HELLO')
    for ind in range(n1,n2):
        #Time intervals are about one minute. So, say an hour interval is 30 on either side
        #print(ind-interval/2)
        Bx = ACE['Bx'][ind-interval/2:ind+interval/2]
        By = ACE['By'][ind-interval/2:ind+interval/2]
        Bz = ACE['Bz'][ind-interval/2:ind+interval/2]

        #calc your n
        
        B = np.transpose(np.array([Bx, By, Bz]))

        M = np.zeros([3,3])
        P = np.zeros([3,3])
        
        B_av = np.array([np.nanmean(B[:,0]), np.nanmean(B[:,1]),np.nanmean(B[:,2])])
        e = B_av/ np.sqrt(B_av[0]**2+B_av[1]**2+B_av[2]**2)
    
        #Create a covariance matrix
        
        for i in range(3):
            for j in range(3):
                M[i,j] = np.nanmean(B[:,i]*B[:,j]) -np.nanmean(B[:,i])*np.nanmean(B[:,j])
                if i == j:
                    P[i,j] = 1 - (e[i]*e[j])
                else:
                    P[i,j] = -1* (e[i]*e[j])
                        
        #Get eigenvalues and eigenvectors
        M = np.dot(np.dot(P,M),P)
        eigenvalues, eigenvectors = LA.eig(M)
        args = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[args]
        eigenvectors = eigenvectors[:,args]
        
        #The vector corresponding to the middle (absolute value) eigenvalue is the minimum variance direction 
        
        front_normal = eigenvectors[:,1]
        
        #The x component of the vector should point towards the sun (positive)
        if front_normal[0] < 0:
            front_normal = -1*front_normal
        
        # Do a test. For the result to be valid, the second smallest eigenvalue should be x (5 for now) times larger than the smallest
        
        
        ratio = eigenvalues[2] / eigenvalues[1]   
    
#        if ratio < 5:
#            print 'MVAB0 ratio is ', ratio
        normals[ind-n1] = front_normal
        ratios[ind-n1] = ratio
        
    vx = ACE['vx'][n1:n2]
    vy = ACE['vy'][n1:n2]
    vz = ACE['vz'][n1:n2]
                    
    shift = np.zeros(len(vx))
    dpos = np.nanmean(ACE['pos'][np.logical_and(ACE['t'] > t1,ACE['t'] < t2)], axis = 0) - (Re* np.nanmean(WIND['pos'][np.logical_and(WIND['t'] > t1,WIND['t'] < t2)], axis = 0))
    for i in range(len(shift)):
        shift[i] = np.dot(normals[i],dpos)/np.dot(normals[i], np.array([vx[i],vy[i],vz[i]]))/60./60./24.
        
    import copy
    t_shifted = copy.deepcopy(ACE['t'][n1:n2])#Overkill? Probably
    t_shifted = t_shifted - shift
    
    p =  ACE[parameter][np.logical_and(ACE['t'] > t1,ACE['t'] < t2)]
    
    p = p[np.argsort(t_shifted)]
    t_shifted = np.sort(t_shifted)
    
    p = p[np.isfinite(t_shifted)]        
    t_shifted = t_shifted[np.isfinite(t_shifted)]   
    #print n
    #print 'MVAB0', shift[0]*24
    #print ''
    return t_shifted, p
#    return t_shifted, p


def getMVAB0Lag(ACE, WIND, t1,t2, interval = 60):
    Re = 6371.
    
    n1 = np.where(ACE['t'] > t1)[0][0]
    n2 = np.where(ACE['t'] < t2)[0][-1]
    
    #Make some variables
    
    normals = np.zeros([len(range(n1,n2)),3])
    ratios = np.zeros([len(range(n1,n2))])
    
    for ind in range(n1,n2):
        #Time intervals are about one minute. So, say an hour interval is 30 on either side
        Bx = ACE['Bx'][ind-interval/2:ind+interval/2]
        By = ACE['By'][ind-interval/2:ind+interval/2]
        Bz = ACE['Bz'][ind-interval/2:ind+interval/2]

        #calc your n
        
        B = np.transpose(np.array([Bx, By, Bz]))

        M = np.zeros([3,3])
        P = np.zeros([3,3])
        
        B_av = np.array([np.nanmean(B[:,0]), np.nanmean(B[:,1]),np.nanmean(B[:,2])])
        e = B_av/ np.sqrt(B_av[0]**2+B_av[1]**2+B_av[2]**2)
    
        #Create a covariance matrix
        
        for i in range(3):
            for j in range(3):
                M[i,j] = np.nanmean(B[:,i]*B[:,j]) -np.nanmean(B[:,i])*np.nanmean(B[:,j])
                if i == j:
                    P[i,j] = 1 - (e[i]*e[j])
                else:
                    P[i,j] = -1* (e[i]*e[j])
                        
        #Get eigenvalues and eigenvectors
        M = np.dot(np.dot(P,M),P)
        eigenvalues, eigenvectors = LA.eig(M)
        args = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[args]
        eigenvectors = eigenvectors[:,args]
        
        #The vector corresponding to the middle (absolute value) eigenvalue is the minimum variance direction 
        
        front_normal = eigenvectors[:,1]
        
        #The x component of the vector should point towards the sun (positive)
        if front_normal[0] < 0:
            front_normal = -1*front_normal
        
        # Do a test. For the result to be valid, the second smallest eigenvalue should be x (5 for now) times larger than the smallest
        
        
        ratio = eigenvalues[2] / eigenvalues[1]   
    
#        if ratio < 5:
#            print 'MVAB0 ratio is ', ratio
        normals[ind-n1] = front_normal
        ratios[ind-n1] = ratio
        
    vx = ACE['vx'][n1:n2]
    vy = ACE['vy'][n1:n2]
    vz = ACE['vz'][n1:n2]
                    
    shift = np.zeros(len(vx))
    flatshift = np.zeros(len(vx))
    dpos = np.nanmean(ACE['pos'][np.logical_and(ACE['t'] > t1,ACE['t'] < t2)], axis = 0) - (Re* np.nanmean(WIND['pos'][np.logical_and(WIND['t'] > t1,WIND['t'] < t2)], axis = 0))
    
    for i in range(len(shift)):
        shift[i] = np.dot(normals[i],dpos)/np.dot(normals[i], np.array([vx[i],vy[i],vz[i]]))/60.
    for i in range(len(shift)):
        flatshift[i] = dpos[0]/vx[i]/60.

    return shift - flatshift
#    return t_shifted, p




def getMVAB0NormalPoints(ACE, WIND, t1,t2, interval = 60):
    Re = 6371.
    
    n1 = np.where(ACE['t'] > t1)[0][0]
    n2 = np.where(ACE['t'] < t2)[0][-1]
    
    #Make some variables
    
    normals = np.zeros([len(range(n1,n2)),3])
    ratios = np.zeros([len(range(n1,n2))])
    
    for ind in range(n1,n2):
        #Time intervals are about one minute. So, say an hour interval is 30 on either side
        Bx = ACE['Bx'][ind-interval/2:ind+interval/2]
        By = ACE['By'][ind-interval/2:ind+interval/2]
        Bz = ACE['Bz'][ind-interval/2:ind+interval/2]

        #calc your n
        
        B = np.transpose(np.array([Bx, By, Bz]))

        M = np.zeros([3,3])
        P = np.zeros([3,3])
        
        B_av = np.array([np.nanmean(B[:,0]), np.nanmean(B[:,1]),np.nanmean(B[:,2])])
        e = B_av/ np.sqrt(B_av[0]**2+B_av[1]**2+B_av[2]**2)
    
        #Create a covariance matrix
        
        for i in range(3):
            for j in range(3):
                M[i,j] = np.nanmean(B[:,i]*B[:,j]) -np.nanmean(B[:,i])*np.nanmean(B[:,j])
                if i == j:
                    P[i,j] = 1 - (e[i]*e[j])
                else:
                    P[i,j] = -1* (e[i]*e[j])
                        
        #Get eigenvalues and eigenvectors
        M = np.dot(np.dot(P,M),P)
        eigenvalues, eigenvectors = LA.eig(M)
        args = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[args]
        eigenvectors = eigenvectors[:,args]
        
        #The vector corresponding to the middle (absolute value) eigenvalue is the minimum variance direction 
        
        front_normal = eigenvectors[:,1]
        
        #The x component of the vector should point towards the sun (positive)
        if front_normal[0] < 0:
            front_normal = -1*front_normal
        
        # Do a test. For the result to be valid, the second smallest eigenvalue should be x (5 for now) times larger than the smallest
        
        
        ratio = eigenvalues[2] / eigenvalues[1]   
    
#        if ratio < 5:
#            print 'MVAB0 ratio is ', ratio
        normals[ind-n1] = front_normal
        ratios[ind-n1] = ratio
        
    return normals, ratios


def getMVAB0NormalPointsWIND(WIND_B, t1,t2, interval = 60):
    Re = 6371.
    
    n1 = np.where(WIND_B['t'] > t1)[0][0]
    n2 = np.where(WIND_B['t'] < t2)[0][-1]
    
    #Make some variables
    
    normals = np.zeros([len(range(n1,n2)),3])
    ratios = np.zeros([len(range(n1,n2))])
    
    for ind in range(n1,n2):
        #Time intervals are about one minute. So, say an hour interval is 30 on either side
        Bx = WIND_B['BGSE'][ind-interval/2:ind+interval/2,0]
        By = WIND_B['BGSE'][ind-interval/2:ind+interval/2,1]
        Bz = WIND_B['BGSE'][ind-interval/2:ind+interval/2,2]

        #calc your n
        
        B = np.transpose(np.array([Bx, By, Bz]))

        M = np.zeros([3,3])
        P = np.zeros([3,3])
        
        B_av = np.array([np.nanmean(B[:,0]), np.nanmean(B[:,1]),np.nanmean(B[:,2])])
        e = B_av/ np.sqrt(B_av[0]**2+B_av[1]**2+B_av[2]**2)
    
        #Create a covariance matrix
        
        for i in range(3):
            for j in range(3):
                M[i,j] = np.nanmean(B[:,i]*B[:,j]) -np.nanmean(B[:,i])*np.nanmean(B[:,j])
                if i == j:
                    P[i,j] = 1 - (e[i]*e[j])
                else:
                    P[i,j] = -1* (e[i]*e[j])
                        
        #Get eigenvalues and eigenvectors
        M = np.dot(np.dot(P,M),P)
        eigenvalues, eigenvectors = LA.eig(M)
        args = np.argsort(eigenvalues)
        eigenvalues = eigenvalues[args]
        eigenvectors = eigenvectors[:,args]
        
        #The vector corresponding to the middle (absolute value) eigenvalue is the minimum variance direction 
        
        front_normal = eigenvectors[:,1]
        
        #The x component of the vector should point towards the sun (positive)
        if front_normal[0] < 0:
            front_normal = -1*front_normal
        
        # Do a test. For the result to be valid, the second smallest eigenvalue should be x (5 for now) times larger than the smallest
        
        
        ratio = eigenvalues[2] / eigenvalues[1]   
    
#        if ratio < 5:
#            print 'MVAB0 ratio is ', ratio
        normals[ind-n1] = front_normal
        ratios[ind-n1] = ratio
        
    return normals, ratios

# Size of basic types (in bytes)
size_logical = 4
size_int = 4
size_double = 8

# For un-aligned data, use '=' (for aligned data set to '')
align = '='

def get_args():
    # Get and parse the command line arguments
    parser = argparse.ArgumentParser(
        description='''Plot .dat files of 1D MPI-AMRVAC simulations''',
        epilog='''Example:
        ./plot_1d.py ...''')
    parser.add_argument('dat_files', nargs='+',
                        type=argparse.FileType('rb'),
                        help='MPI-AMRVAC .dat files to read')
    return parser.parse_args()

def read_header(dat):
    h = {}

    fmt = align + 'iiiiiiiiiid'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['version'], h['offset_tree'], h['offset_blocks'], h['nw'],
     h['ndir'], h['ndim'], h['levmax'], h['nleafs'], h['nparents'],
     h['it'], h['time']] = hdr

    fmt = align + 2 * 'd' + 2 * 'i'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    [h['xmin'], h['xmax'], h['domain_nx'], h['block_nx']] = hdr

    # Read w_names
    w_names = []
    for i in range(h['nw']):
        fmt = align + 16 * 'c'
        hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w_names.append(''.join(hdr).strip())
    h['w_names'] = w_names

    fmt = align + 16 * 'c'
    hdr = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    h['physics_type'] = ''.join(hdr).strip()

    fmt = align + 'i'
    n_pars, = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    fmt = align + n_pars * 'd'
    vals = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    fmt = align + n_pars * 16 * 'c'
    names = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
    # Split and join the name strings (from one character array)
    names = [''.join(names[i:i+16]).strip() for i in range(0, len(names), 16)]

    for val, name in zip(vals, names):
        h[name] = val

    return h

def get_data_1d(dat):
    h = read_header(dat)
    nw = h['nw']
    nx = h['block_nx']
    nleafs = h['nleafs']
    nparents = h['nparents']

    # Read tree info. Skip 'leaf' array
    dat.seek(h['offset_tree'] + (nleafs+nparents) * size_logical)

    # Read block levels
    fmt = align + nleafs * 'i'
    block_lvls = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Read block indices
    fmt = align + nleafs * 'i'
    block_ixs = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

    # Determine coarse grid spacing
    dx0 = (h['xmax'] - h['xmin']) / h['domain_nx']

    ### Start reading data blocks
    dat.seek(h['offset_blocks'])

    out_data = np.zeros((nx * nleafs, nw+1))

    for i in range(nleafs):
        lvl = block_lvls[i]
        ix = block_ixs[i]

        x0 = (ix-1) * nx * dx0 * 0.5**(lvl-1) + dx0 * 0.5**lvl + h['xmin']
        x1 = ix * nx * dx0 * 0.5**(lvl-1) - dx0 * 0.5**lvl + h['xmin']
        x = np.linspace(x0, x1, nx)

        # Read number of ghost cells
        fmt = align + 2 * 'i'
        [gc_lo, gc_hi] = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))

        # Read actual data
        block_size = gc_lo + nx + gc_hi
        fmt = align + block_size * nw * 'd'
        d = struct.unpack(fmt, dat.read(struct.calcsize(fmt)))
        w = np.reshape(d, [block_size, nw], 'F') # Fortran ordering
        out_data[i*nx:(i+1)*nx, 0] = x
        out_data[i*nx:(i+1)*nx, 1:] = w[gc_lo:gc_lo+nx, :]

    return h, out_data

def get_primitive_data(cc, wnames, h):
    pp = np.copy(cc)
    prim_names = wnames[:]
    ndir = h['ndir']
    i_rho = wnames.index('rho')

    # Convert momentum
    try:
        i_m1 = wnames.index('m1')
        for i in range(ndir):
            pp[:, i_m1+i] = pp[:, i_m1+i] / pp[:, i_rho]
            prim_names[i_m1+i] = 'v' + str(i+1)
    except ValueError:
        pass

    # Convert energy
    try:
        i_e = wnames.index('e')
        prim_names[i_e] = 'p'
        kin_en = 0.5 * pp[:, i_rho] * np.sum(pp[:, i_m1:i_m1+ndir]**2, axis=1)

        if h['physics_type'] == 'hd':
            pp[:, i_e] = (h['gamma'] - 1.0) * (pp[:, i_e] - kin_en)
        elif h['physics_type'] == 'mhd':
            i_b1 = wnames.index('b1')
            mag_en = 0.5 * np.sum(pp[:, i_b1:i_b1+ndir]**2, axis=1)
            pp[:, i_e] = (h['gamma'] - 1.0) * (pp[:, i_e] - kin_en - mag_en)
        else:
            print("Unknown physics type, cannot convert to pressure")
    except ValueError:
        pass

    return [pp, prim_names]

def read_files(fpath):

    import os
    #fpath = 'C:/Users/Taylor/VMshare/sw_sim/'
    names = os.listdir(fpath)
    files = []
    for name in names:
        if name[-3:] == 'dat':
            files.append(fpath+name)

    if files == []:
        print('No data files :(')
        return 0,0,0

    #sort names
    files = np.sort(files)
    
    all_cons = []
    all_prim = []
    all_times = []
    
    # Read in all data
    for f in files:
        h, data_1d = get_data_1d(open(f, 'rb'))
        all_cons.append(data_1d)
        cons_names = ['x'] + h['w_names']
        prim_1d, prim_names = get_primitive_data(data_1d, cons_names, h)
        all_prim.append(prim_1d)
        all_times.append(h['time']) 
    qlength_unit = 10**3
    qmass_unit = 10**-12 
    #Apply conversion
    for i in range(len(all_times)):
        all_prim[i][:,1] = all_prim[i][:,1]/ 1.e6 / (1.6726e-27/qmass_unit) / qlength_unit**3 #/1.6726
        all_prim[i][:,6] = all_prim[i][:,6]/1.e-9 * np.sqrt(1.257e-6 * qmass_unit/qlength_unit) # /28.209
        all_prim[i][:,7] = all_prim[i][:,7]/1.e-9 * np.sqrt(1.257e-6 * qmass_unit/qlength_unit) # / 28.209
        all_prim[i][:,8] = all_prim[i][:,8]/1.e-9 * np.sqrt(1.257e-6 * qmass_unit/qlength_unit) # / 28.209
        all_prim[i][:,5] = all_prim[i][:,5] #/ all_prim[i][:,1] *1.e6 / qlength_unit**3 / 1.38065e-23 * (qmass_unit * qlength_unit**2)  #/(all_prim[i][:,1]* 1.6726) / 0.008254 
        all_prim[i][:,2] = all_prim[i][:,2]*1.e3 / qlength_unit
        all_prim[i][:,3] = all_prim[i][:,3]*1.e3 / qlength_unit
        all_prim[i][:,4] = all_prim[i][:,4]*1.e3 / qlength_unit

    
    return all_times, prim_names,all_prim


def read_files_cons(fpath):

    import os
    #fpath = 'C:/Users/Taylor/VMshare/sw_sim/'
    names = os.listdir(fpath)
    files = []
    for name in names:
        if name[-3:] == 'dat':
            files.append(fpath+name)

    if files == []:
        print('No data files :(')
        return 0,0,0

	#sort names
    files = np.sort(files)
		
    all_cons = []
    all_times = []
    
    # Read in all data
    for f in files:
        h, data_1d = get_data_1d(open(f, 'rb'))
        all_cons.append(data_1d)
        cons_names = ['x'] + h['w_names']
        all_times.append(h['time']) 

    
    return all_times, cons_names,all_cons