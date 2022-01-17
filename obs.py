# -*- coding: utf-8 -*-
"""
1. Create obs with random error [Save: obs.txt]
2. Test the error covariance

"""
import numpy as np
import matplotlib.pyplot as plt
import os

path = 'D:\\Desktop file\\All Files\\NTU課程資料\\大氣系課程\\資料同化-連國淵老師\\Project\\code'
os.chdir(path)

# read nature run (x_t)
x_t = np.genfromtxt("txt_file_new/x_t.txt")

# Numbers of timesteps
nT = 400

#%% create obs error 
mu = 0
sigma = 0.1
e1 = np.random.normal(mu, sigma, (nT,40))     # original obs, higher obs error
e2 = np.random.normal(mu, 0.5*sigma, (nT,40)) # newer obs, lower obs error

# error for testing (use e1 for example)
err1 = e1[:,0]
err2 = e1[:,2]

# Save files
np.savetxt('txt_file_new/e1.txt', e1)
np.savetxt('txt_file_new/e2.txt', e2)

#%% Extract specific data by obs arrangement (Exp 1.)

# Combine obs and error [size=(400,40)]
obs_1 = x_t[0:nT,:] + e1  # Use old obs for all time steps.
obs_2 = x_t[0:nT,:] + e2  # Use new obs for all time steps.
obs_3 = np.vstack((obs_1[0:50,:],obs_2[50:nT,:])) # 50 old, 350 new.
'''different ratio could apply to obs_3.'''

# 20 evenly distributed stations (0,2,4,6,...,38)
L = np.int(20)
# obs (400, 20)
obs1 = obs_1[:,::2]
obs2 = obs_2[:,::2]
obs3 = obs_3[:,::2]
# error covariance metrix R (20,20)
R1 = np.zeros((L,L))  # old obs R
R2 = np.zeros((L,L))  # new obs R
np.fill_diagonal(R1, 0.01)
np.fill_diagonal(R2, 0.0025)

# Observation operator
H = np.zeros((L,40))
i = 0
k = 0
for i in range(20):
    H[i,k] += 1
    k += 2

# Save files
np.savetxt('txt_file_new/exp1_obs1.txt', obs1)
np.savetxt('txt_file_new/exp1_obs2.txt', obs2)
np.savetxt('txt_file_new/exp1_obs3.txt', obs3)
np.savetxt('txt_file_new/exp1_R1.txt',   R1)
np.savetxt('txt_file_new/exp1_R2.txt',   R2)
np.savetxt('txt_file_new/exp1_H.txt',    H)

#%% Extract specific data by obs arrangement (Exp 2.)
'''
Case1: 20 evenly distributed stations (0,2,4,6,...,38)
Case2: When tt>50, 6 places move toward right and left. (4 to 3, 6 to 7, ...)
Case3: For all tt, 6 places move toward right and left. (4 to 3, 6 to 7, ...)
''' 

# Combine obs and error [size=(400,40)]
obs = x_t[0:nT,:] + e1  # Use old obs for all time steps.

L = np.int(20)
# obs (400, 20)
obs1 = obs[:,::2]
obs2 = np.zeros((nT,20))
obs3 = np.zeros((nT,20))
location = [0,2,3,7,8,10,12,14,15,19,20,22,24,26,27,31,32,34,36,38]
for i in range(20):
    obs2[ 0:50,i] += obs[ 0:50,i*2]  # original arrangement
    obs2[50:nT,i] += obs[50:nT,location[i]]
    obs3[    :,i] += obs[    :,location[i]]

# error covariance metrix R (20,20)
R = np.zeros((L,L))  # old obs R
np.fill_diagonal(R, 0.01)

# Observation operator H (20,40)
H1 = np.zeros((L,40))
H2 = np.zeros((L,40)) # for case 2 and 3
i = 0
k = 0
for i in range(20):
    H1[i,k] += 1
    H2[i,location[i]] += 1
    k += 2

# Save files
np.savetxt('txt_file_new/exp2_obs1.txt', obs1)
np.savetxt('txt_file_new/exp2_obs2.txt', obs2)
np.savetxt('txt_file_new/exp2_obs3.txt', obs3)
np.savetxt('txt_file_new/exp2_R.txt',    R)
np.savetxt('txt_file_new/exp2_H1.txt',   H1)
np.savetxt('txt_file_new/exp2_H2.txt',   H2)

#%% Extract specific data by obs arrangement (Exp 3.)
'''
Case1: Evenly distributed pairs of stations (0,1,4,5,8,9,...,32,33,36,37)
Case2: Moving right hand side at every time step.
''' 
e1 = np.genfromtxt("txt_file_new/e1.txt")

# Combine obs and error [size=(400,40)]
obs = x_t[0:nT,:] + e1  # Use old obs for all time steps.

L = np.int(20)
# obs (400, 20)
loc1_exp3 = [0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33,36,37]
loc2_exp3 = [1,2,5,6,9,10,13,14,17,18,21,22,25,26,29,30,33,34,37,38]
loc3_exp3 = [2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35,38,39]
loc4_exp3 = [0,3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32,35,36,39]
obs1 = np.zeros((nT,20))
obs2 = np.zeros((nT,20))
i = 0
j = 0
k = 0
for i in range(20):
    obs1[:,i] += obs[:,loc1_exp3[i]]
for j in range(nT):
    if j%4 == 0:
        for k in range(20):
            obs2[j,k] = obs[j,loc1_exp3[k]]
    elif j%4 == 1:
        for k in range(20):
            obs2[j,k] = obs[j,loc2_exp3[k]]
    elif j%4 == 2:
        for k in range(20):
            obs2[j,k] = obs[j,loc3_exp3[k]]
    elif j%4 == 3:
        for k in range(20):
            obs2[j,k] = obs[j,loc4_exp3[k]]

# error covariance metrix R (20,20)
R = np.zeros((L,L))  # old obs R
np.fill_diagonal(R, 0.01)

# Observation operator H (20,40) => H will keep changing at each time step.
H1 = np.zeros((L,40))
H2 = np.zeros((L,40))
H3 = np.zeros((L,40))
H4 = np.zeros((L,40))
i = 0
for i in range(20):
    H1[i,loc1_exp3[i]] += 1
    H2[i,loc2_exp3[i]] += 1
    H3[i,loc3_exp3[i]] += 1
    H4[i,loc4_exp3[i]] += 1

# Save files
np.savetxt('txt_file_new/exp3_obs1.txt', obs1)
np.savetxt('txt_file_new/exp3_obs2.txt', obs2)
np.savetxt('txt_file_new/exp3_R.txt',    R)
np.savetxt('txt_file_new/exp3_H1.txt',   H1)
np.savetxt('txt_file_new/exp3_H2.txt',   H2)
np.savetxt('txt_file_new/exp3_H3.txt',   H3)
np.savetxt('txt_file_new/exp3_H4.txt',   H4)

#%% Testing by scatter plot

def scatter_hist(x, y, ax, ax_histx, ax_histy):
    
    colorname = 'mediumorchid'
    
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y, ec='indigo', c=colorname)

    # now determine nice limits by hand:
    binwidth = 0.015
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins, ec='black', fc=colorname)
    ax_histy.hist(y, bins=bins, orientation='horizontal', ec='black', fc=colorname)
    ax_histx.set_ylim(0,35)
    ax_histy.set_xlim(0,35)
    ax_histx.set_yticks([0,5,10,15,20,25,30,35])
    ax_histy.set_xticks([0,5,10,15,20,25,30,35])


# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.01
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]
# start with a square Figure
fig = plt.figure(figsize=(8, 8))
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

# use the previously defined function
e_obs1 = err1
e_obs2 = err2
scatter_hist(e_obs1, e_obs2, ax, ax_histx, ax_histy)

plt.savefig('plot_new/test_obs_error.png',dpi=300)
plt.tight_layout()
plt.show()



