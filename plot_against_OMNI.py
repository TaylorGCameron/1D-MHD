# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 15:05:08 2018

@author: Taylor
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
import datetime
import amrvac_support as vac
import csv
import os
from ai import cdas



data = np.load('output.npy')

t = data['t']
n = data['n']
vx = data['vx']
Bz = data['Bz']
myFmt = mdate.DateFormatter('%H:%M')

OMNI = cdas.get_data('sp_phys','OMNI_HRO_1MIN', mdate.num2date(t[0]),mdate.num2date(t[-1]),
['BX_GSE','BY_GSE','BZ_GSE', 'Vx','Vy','Vz','proton_density','T','Pressure','E','AE_INDEX','BSN_x','BSN_y','BSN_z']
)

OMNI['VX_VELOCITY,_GSE'][OMNI['VX_VELOCITY,_GSE'] > 9999] = np.nan
OMNI['BZ,_GSE'][OMNI['BZ,_GSE'] > 9999] = np.nan
OMNI['PROTON_DENSITY'][OMNI['PROTON_DENSITY'] > 999] = np.nan

plt.figure(2)
myFmt = mdate.DateFormatter('%H:%M')


ax1 = plt.subplot(311)
ax1.plot(mdate.num2date(t), n)
ax1.plot(OMNI['EPOCH_TIME'], OMNI['PROTON_DENSITY'])
#ax1.set_ylim(0,50)
ax1.set_ylabel('Density')

ax2 = plt.subplot(312)
ax2.plot(mdate.num2date(t), vx)
ax2.plot(OMNI['EPOCH_TIME'], OMNI['VX_VELOCITY,_GSE'])
#ax2.set_ylim(-600,-300)
ax2.set_ylabel('Vx')

ax3 = plt.subplot(313)
ax3.plot(mdate.num2date(t), Bz)
ax3.plot(OMNI['EPOCH_TIME'], OMNI['BZ,_GSE'])
#ax3.set_ylim(-20,20)
ax3.xaxis.set_major_formatter(myFmt)
ax3.set_ylabel('IMF Bz')
