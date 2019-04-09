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

data_prefix = 'data/'


#Load in model data
t, names, data = vac.read_files(os.getcwd()+'/')
if t == 0:
    import sys 
    sys.exit()
    


n = 64

#Make the plot
    
fig1 = plt.figure(3, figsize = (12,12))
ax1 = plt.subplot(411)
ax1.set_ylabel('')
ax2 = plt.subplot(412)
ax2.set_ylabel('rho')
ax3 = plt.subplot(413)
ax3.set_ylabel('B')
ax4 = plt.subplot(414)
ax4.set_ylabel('p')


ax1.plot(data[n][:,0],data[n][:,1])
ax2.plot(data[n][:,0],data[n][:,2])
ax2.plot(data[n][:,0],data[n][:,3])
ax2.plot(data[n][:,0],data[n][:,4])
ax3.plot(data[n][:,0],data[n][:,6])
ax3.plot(data[n][:,0],data[n][:,7])
ax3.plot(data[n][:,0],data[n][:,8])
ax4.plot(data[n][:,0],data[n][:,5])


#Okay, so based on this, there's a problem occurring in the last frame, where the pressure is suddenly dropping, which I assume sets off the crash.