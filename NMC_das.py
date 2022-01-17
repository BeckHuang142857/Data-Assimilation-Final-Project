# -*- coding: utf-8 -*-
"""
das for NMC: create analysis using P_b = R

[Use the same initial condition.]
1. perfect model & e1   => pm1
2. perfect model & e2   => pm2
3. imperfect model & e1 => im1
4. imperfect model & e2 => im2

"""

import os
import numpy as np
from numpy.linalg import inv, multi_dot
from scipy.integrate import ode
import lorenz96
import matplotlib.pyplot as plt
#from setting_perfect_model import *
from setting_imperfect_model import *

path = 'D:\\Desktop file\\All Files\\NTU課程資料\\大氣系課程\\資料同化-連國淵老師\\Project\\code'
os.chdir(path)

#%% Load files

# load initial condition
x_a_init = np.genfromtxt('txt_file_new/x_a_init_perfect.txt')

# load observations ( 2 kinds of obs )
#y_o_save = np.genfromtxt('txt_file_new/NMC_full_obs_e1.txt') 
y_o_save = np.genfromtxt('txt_file_new/NMC_full_obs_e2.txt') 

# Error covariance metrix R (20,20)
R = np.zeros((40,40))  
np.fill_diagonal(R, 0.01)

# Observation operator H
H = np.zeros((40,40))  
np.fill_diagonal(H, 1)

# Background error covariance P_b 
P_b = R

# initial x_b: no values at the initial time (assign NaN)
x_b_save = np.full((1,N), np.nan, dtype='f8')

# initial x_a: from x_a_init
x_a_save = np.array([x_a_init])

#%% Data Assimilation Cycle

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

    # innovation
    y_b = np.dot(H, x_b)
    d = y_o - y_b
    
    # analysis scheme (OI)
    inverse = inv(multi_dot([H, P_b, H.transpose()]) + R)
    K = multi_dot([P_b, H.transpose(), inverse])
    x_a = x_b + np.dot(K, d)

    x_a_save = np.vstack([x_a_save, x_a.transpose()])
    tt += 1

name_list = ['pm1','pm2','im1','im2']
name = name_list[3]

# save background and analysis data
np.savetxt('txt_file_new/NMC_x_b_'+name+'.txt', x_b_save)
np.savetxt('txt_file_new/NMC_x_a_'+name+'.txt', x_a_save)


#%%  Test

x_t = np.genfromtxt("txt_file_new/x_t.txt")
diff = x_a_save - x_t

plt.plot(diff[:,0])
plt.show()

