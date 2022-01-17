# -*- coding: utf-8 -*-
"""
The data assimilation system (for Exp 2)
Load:
  x_a_init.txt
  exp2_obs1.txt
  exp2_obs2.txt
  exp2_obs3.txt
  exp2_R.txt
  exp2_H1.txt
  exp2_H2.txt
Save:
  exp2_?_x_b.txt
  exp2_?_x_a.txt
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

# load observations ( 3 kinds of obs )
''' Change this for each run '''
y_o_save = np.genfromtxt('txt_file_new/exp2_obs1.txt') # 20 evenly distributed stations (0,2,4,6,...,38)
#y_o_save = np.genfromtxt('txt_file_new/exp2_obs2.txt') # When tt>50, 6 places move toward right and left. (4 to 3, 6 to 7, ...)
#y_o_save = np.genfromtxt('txt_file_new/exp2_obs3.txt') # For all tt, 6 places move toward right and left. (4 to 3, 6 to 7, ...)

# Error covariance metrix R (20,20)
R = np.genfromtxt('txt_file_new/exp2_R.txt')

# Observation operator H
H1 = np.genfromtxt('txt_file_new/exp2_H1.txt') # original arrangement
H2 = np.genfromtxt('txt_file_new/exp2_H2.txt') # new arrangement
H  = H1

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
exp = 1

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
    
    # Change H when tt > 50 in exp2_2
    if exp==2 and tt>50:
        H = H2
    elif exp==3:
        H = H2
    
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
np.savetxt('txt_file_new/exp2_'+str(exp)+'_x_b_'+mn+'.txt', x_b_save)
np.savetxt('txt_file_new/exp2_'+str(exp)+'_x_a_'+mn+'.txt', x_a_save)
print('Exp 2-',exp,' ',mn,' Finish')









