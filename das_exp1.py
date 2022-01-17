# -*- coding: utf-8 -*-
"""
The data assimilation system (for Exp 1)
Load:
  x_a_init.txt
  exp1_obs1.txt
  exp1_obs2.txt
  exp1_obs3.txt
  exp1_R1.txt
  exp1_R2.txt
  exp1_H.txt
Save:
  exp1_?_x_b.txt
  exp1_?_x_a.txt
"""
import os
import numpy as np
from numpy.linalg import inv, multi_dot
from scipy.integrate import ode
import lorenz96

pm = False  # Using "Perfect model" or "Imperfect model"
if pm==True:
    from setting_perfect_model import *
else:
    from setting_imperfect_model import *

path = 'D:\\Desktop file\\All Files\\NTU課程資料\\大氣系課程\\資料同化-連國淵老師\\Project\\code'
os.chdir(path)

#%% Load files

# load initial condition
x_a_init = np.genfromtxt('txt_file_new/x_a_init_perfect.txt')

# load observations ( 3 kinds of obs )
''' Change this for each run '''
#y_o_save = np.genfromtxt('txt_file_new/exp1_obs1.txt') # Use old obs for all time steps.
#y_o_save = np.genfromtxt('txt_file_new/exp1_obs2.txt') # Use new obs for all time steps.
y_o_save = np.genfromtxt('txt_file_new/exp1_obs3.txt') # 50 old, 350 new.

# Error covariance metrix R (20,20)
R1 = np.genfromtxt('txt_file_new/exp1_R1.txt') # for old obs
R2 = np.genfromtxt('txt_file_new/exp1_R2.txt') # for new obs

# Observation operator H
H = np.genfromtxt('txt_file_new/exp1_H.txt')

# Background error covariance P_b (NMC method)
if pm==True:
    P_b_1 = np.genfromtxt('txt_file_new/NMC_P_b_pm1.txt')
    P_b_2 = np.genfromtxt('txt_file_new/NMC_P_b_pm2.txt')
else:
    P_b_1 = np.genfromtxt('txt_file_new/NMC_P_b_im1.txt')
    P_b_2 = np.genfromtxt('txt_file_new/NMC_P_b_im2.txt')

# initial x_b: no values at the initial time (assign NaN)
x_b_save = np.full((1,N), np.nan, dtype='f8')

# initial x_a: from x_a_init
x_a_save = np.array([x_a_init])

#%% Data Assimilation Cycle

''' Change this index for each run '''
exp = 4
# Choose R for each run
if exp==1:   # Use old obs for all time steps.
    R = R1 
    P_b = P_b_1
elif exp==2: # Use new obs for all time steps.
    R = R2
    P_b = P_b_2
elif exp==3: # 50 old, 350 new. (use correct error in model)
    R = R1
    P_b = P_b_1
elif exp==4: # 50 old, 350 new. (still use old error in model)
    R = R1
    P_b = P_b_1

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
    
    # Change R and P_b when tt > 50 in exp1_3
    if exp==3 and tt>50:
       R = R2
       P_b = P_b_2 
    
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
np.savetxt('txt_file_new/exp1_'+str(exp)+'_x_b_'+mn+'.txt', x_b_save)
np.savetxt('txt_file_new/exp1_'+str(exp)+'_x_a_'+mn+'.txt', x_a_save)
print('Exp 1-',exp,' Finish')







