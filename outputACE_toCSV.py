# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 12:48:48 2018

@author: Taylor
"""
import numpy as np
import datetime
import matplotlib.dates as mdate
import amrvac_support as vac
import csv
import matplotlib.pyplot as plt


######################################################### 
#Read in model_config.par
#########################################################
lines = []
f = open('model_config.par')
pars= f.read()
f.close()

t1 = np.array(vac.get_keyword(pars,'interval_start').split(',')).astype('int')
t2 = np.array(vac.get_keyword(pars,'interval_end').split(',')).astype('int')

t1_mvab0 = np.array(vac.get_keyword(pars,'mvab0_interval_start').split(',')).astype('int')
t2_mvab0 = np.array(vac.get_keyword(pars,'mvab0_interval_end').split(',')).astype('int')

time_integrator = vac.get_keyword(pars,'time_integrator')
flux_scheme = vac.get_keyword(pars,'flux_scheme')
limiter = vac.get_keyword(pars,'limiter')
courantpar = vac.get_keyword(pars,'courantpar')

#########################################################
#Get start and end times
#########################################################
year = t1[0]
t1 = mdate.date2num(datetime.datetime(t1[0],t1[1],t1[2],t1[3],t1[4],t1[5]) - datetime.timedelta(0./24.))
t2 = mdate.date2num(datetime.datetime(t2[0],t2[1],t2[2],t2[3],t2[4],t2[5]))

t1_mvab0 = mdate.date2num(datetime.datetime(t1_mvab0[0], t1_mvab0[1], t1_mvab0[2], t1_mvab0[3], t1_mvab0[4],t1_mvab0[5]))
t2_mvab0 = mdate.date2num(datetime.datetime(t2_mvab0[0], t2_mvab0[1], t2_mvab0[2], t2_mvab0[3], t2_mvab0[4],t2_mvab0[5]))

#########################################################
#Read in, rotate and process ACE data
#########################################################

swe_data, mfi_data = vac.pull_ACE(t1, t2)
mfi_data = vac.average_mfi(swe_data, mfi_data)
ACE = vac.combine_ACE(swe_data, mfi_data)
ACE = vac.clean_ACE(ACE)

#Get the MVAB-0 normal
n = np.array([1,0,0])
normal, ratio = vac.getMVAB0Normal(t1_mvab0, t2_mvab0)

#Rotate ACE data into coordinates aligned with this normal
axis = vac.get_axis(np.array([1,0,0]), normal)
angle = vac.get_angle(np.array([1,0,0]), normal)
ACE =  vac.rotate_ACE(ACE, -1.*angle , axis)



#########################################################
#Write to files
#########################################################

#Write the data to a csv to be read in by AMRVAC
with open('data.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    for i in range(len(ACE['t'])):
        writer.writerow([(ACE['t'][i] - t1)*24*60.*60., ACE['n'][i], ACE['vx'][i], ACE['vy'][i], ACE['vz'][i], ACE['Bx'][i], ACE['By'][i], ACE['Bz'][i], ACE['T'][i]])

#Calculate model length as 1.15 times the distance between ACE and the Earth.
model_length = round(np.max(ACE['pos'][:,0])*1.15, -4)
interval_length = int((t2-t1)*24.*60.*60.) #hours

#Write in the interval, MVAB-0 interval, the phase front normal and model length to file for use by other scripts.
with open('interval_information.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow([t1,t2])
    writer.writerow([t1_mvab0,t2_mvab0])
    writer.writerow(normal)
    writer.writerow([model_length])


#Modify amrvac parameter file
values = [
             ['xprobmax1',('%.2E' % model_length).replace('E+','d')],
             ['time_max',str(interval_length)+'d0'],
             ['time_integrator',time_integrator],
             ['flux_scheme',flux_scheme],
             ['limiter',limiter],
             ['courantpar',courantpar]
         ]

vac.write_keywords(values)

