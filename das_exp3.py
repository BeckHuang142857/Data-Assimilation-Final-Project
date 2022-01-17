# -*- coding: utf-8 -*-
"""
The data assimilation system (for Exp 3)
Load:
  x_a_init.txt
  exp3_obs1.txt
  exp3_obs2.txt
  exp3_R.txt
  exp3_H1.txt
  exp3_H2.txt
  exp3_H3.txt
  exp3_H4.txt
Save:
  exp3_?_x_b.txt
  exp3_?_x_a.txt
"""
import os
import numpy as np
from numpy.linalg import inv, multi_dot
from scipy.integrate import ode
import lorenz96

pm = False  # Using "Perfect model"(pm) or "Imperfect model"(im)
if pm==True:
    from setting_perfect_model import *
else:
    from setting_imperfect_model import *

path = 'D:\\Desktop file\\All Files\\NTU課程資料\\大氣系課程\\資料同化-連國淵老師\\Project\\code'
os.chdir(path)

#%% Load files

# load initial condition
x_a_init = np.genfromtxt('txt_file_new/x_a_init_perfect.txt')

# load observations ( 2 kinds of obs )
''' Change this for each run '''
#y_o_save = np.genfromtxt('txt_file_new/exp3_obs1.txt') # Use static obs for all time steps.
y_o_save = np.genfromtxt('txt_file_new/exp3_obs2.txt') # Use moving obs for all time steps.

# Error covariance metrix R (20,20)
R = np.genfromtxt('txt_file_new/exp3_R.txt') 

# Observation operator H
H1 = np.genfromtxt('txt_file_new/exp3_H1.txt')  
H2 = np.genfromtxt('txt_file_new/exp3_H2.txt')
H3 = np.genfromtxt('txt_file_new/exp3_H3.txt')
H4 = np.genfromtxt('txt_file_new/exp3_H4.txt')
H  = H1                                    # default

# Background error covariance P_b (NMC method)
if pm==True:
    P_b = np.genfromtxt('txt_file_new/NMC_P_b_pm1.txt')
else:
    P_b = np.genfromtxt('txt_file_new/NMC_P_b_im1.txt')

# initial x_b: no values at the initial time (assign NaN)
x_b_save = np.full((1,N), np.nan, dtype='f8')

# initial x_a: from x_a_init
x_a_save = np.array([x_a_init])

#%% Data Assimilation Cycle

''' Change this index for each run '''
exp = 2

tt = 1
while tt <= nT:
    tts = tt - 1
    Ts = tts * dT  # forecast start time
    Ta = tt  * dT  # forecast end time (DA analysis time)
    print('Cycle =', tt, ', Ts =', round(Ts, 10), ', Ta =', round(Ta, 10))

    #--------------
    # forecast step
    #--------------

    solver = ode(lorenz96.f).set_integrator('dopri5')
    solver.set_initial_value(x_a_save[tts], Ts).set_f_params(F)
    solver.integrate(Ta)
    x_b_save = np.vstack([x_b_save, [solver.y]])

    #--------------
    # analysis step
    #--------------

    # background
    x_b = x_b_save[tt].transpose()

    # observation
    y_o = y_o_save[tt-1,:]
    
    # Change H for each tt
    if exp==2:
        if tt%4 == 1:
            H = H1
        elif tt%4==2:
            H = H2
        elif tt%4==3:
            H = H3
        elif tt%4==0:
            H = H4
    
    # innovation
    y_b = np.dot(H, x_b)
    d = y_o - y_b  
    
    # analysis scheme (OI)
    inverse = inv(multi_dot([H, P_b, H.transpose()]) + R)
    K = multi_dot([P_b, H.transpose(), inverse])
    x_a = x_b + np.dot(K, d)

    x_a_save = np.vstack([x_a_save, x_a.transpose()])
    tt += 1

# save background and analysis data
if pm==True:
    mn = 'pm' # mn = model name
else:
    mn = 'im'
np.savetxt('txt_file_new/exp3_'+str(exp)+'_x_b_'+mn+'.txt', x_b_save)
np.savetxt('txt_file_new/exp3_'+str(exp)+'_x_a_'+mn+'.txt', x_a_save)
print('Exp 3-',exp,' ',mn,' Finish')




