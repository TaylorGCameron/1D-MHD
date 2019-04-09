# -*- coding: utf-8 -*-
"""
Created on Wed Aug 08 15:00:57 2018

@author: Taylor
"""

#plot_amrvac_againt_WIND

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import datetime
import amrvac_support as vac
import csv
import os



#########################################################
#load in data
#########################################################

#Load in data from CSVs
extra = []
with open('interval_information.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        extra.append(row)   
t0 = float(extra[0][0]) #beginning of input data, also simulation
t1 = float(extra[0][1]) #end of input data
t0_mvab0 = float(extra[1][0]) #beginning of mvab0 interval
t1_mvab0 = float(extra[1][1]) #end of mvab0 interval
normal = np.array(extra[2]).astype('float')
model_length = float(extra[3][0])

axis = vac.get_axis(np.array([1,0,0]), normal)
angle = vac.get_angle(np.array([1,0,0]), normal)


#Load in ACE data
swe_data, mfi_data = vac.pull_ACE(t0, t1)
mfi_data = vac.average_mfi(swe_data, mfi_data)
ACE = vac.combine_ACE(swe_data, mfi_data)
ACE = vac.clean_ACE(ACE)


#load in bowshock nose location from OMNI
BSN = vac.pull_OMNI_BSN(t0,t1)

#Load in model data

t, names, data = vac.read_files(os.getcwd()+'/')

#########################################################
#Process the data
#########################################################

x = np.zeros(len(t))
n = np.zeros(len(t))
vx = np.zeros(len(t))
vy = np.zeros(len(t))
vz = np.zeros(len(t))

Bx = np.zeros(len(t))
By = np.zeros(len(t))
Bz = np.zeros(len(t))

na = np.zeros(len(t))
vxa = np.zeros(len(t))
Bza = np.zeros(len(t))


vwx = np.zeros(len(t))
ACEx= np.zeros(len(t))
WINDx= np.zeros(len(t))



#generate virtual BSN time series
t = np.array(t)/60./60./24. + t0 #turn model t into real time t

Apos = ACE['pos']
bsnx = BSN['bsn_x_gse'] * 6371.
bsny = BSN['bsn_y_gse'] * 6371.
bsnz = BSN['bsn_z_gse'] * 6371.
BSNpos = np.transpose([bsnx,bsny,bsnz])
for i in range(len(Apos)):
    Apos[i] = vac.rotate(Apos[i], -1*angle, axis)
for i in range(len(BSNpos)):
    BSNpos[i] = vac.rotate(BSNpos[i], -1*angle, axis)
Apos_t = ACE['t']
BSNpos_t = BSN['t']



#Interpolate and sample model data to get virtual time series
for i in range(len((t))):
    dist_x = np.interp(t[i], Apos_t, Apos[:,0]) - np.interp(t[i], BSNpos_t, BSNpos[:,0])
    #dist_x = 1133197.4791891947 #avg dist_x in shock coords
    #model_length = 1510000. #This is the end of the grid
    virtual_BSN_x = model_length - dist_x
    vwx[i] = virtual_BSN_x
    n[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,1])
    vx[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,2])
    vy[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,3])
    vz[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,4])

    Bx[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,6])
    By[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,7])
    Bz[i] = np.interp(virtual_BSN_x, data[i][:,0], data[i][:,8])

    na[i] = data[i][-1,1]
    vxa[i] = data[i][-1,2]
    Bza[i] = data[i][-1,8]



#Ratate v and B back to GSE
    
v = np.transpose([vx, vy, vz])
vr = np.empty([len(v), 3])
for i in range(len(v)):
    vr[i] = vac.rotate(v[i], angle, axis)  
vx = vr[:,0]
vy = vr[:,1]
vz = vr[:,2]

B = np.transpose([Bx, By, Bz])
Br = np.empty([len(B), 3])
for i in range(len(B)):
    Br[i] = vac.rotate(B[i], angle, axis)  
Bz = Br[:,2]
Bx = Br[:,0]
By = Br[:,1]


#Make a structure to hold the shifted data
dtype = np.dtype([('t', 'f8'),('n', 'f8'),('vx', 'f8'),('vy', 'f8'),('vz', 'f8'),('Bx', 'f8'),('By', 'f8'),('Bz', 'f8')])
shifted_data = np.ndarray(len(t), dtype = dtype)        

shifted_data['t'] = t
shifted_data['n'] = n
shifted_data['vx'] = vx
shifted_data['vy'] = vy
shifted_data['vz'] = vz
shifted_data['Bx'] = Bx
shifted_data['By'] = By
shifted_data['Bz'] = Bz


#save the data
np.save('output.npy', shifted_data)



