# -*- coding: utf-8 -*-
"""
Plot results for Exp.1
"""
import numpy as np
import matplotlib.pyplot as plt
import os

path = 'D:\\Desktop file\\All Files\\NTU課程資料\\大氣系課程\\資料同化-連國淵老師\\Project\\code'
os.chdir(path)

# read nature run (x_t)
x_t = np.genfromtxt("txt_file_new/x_t.txt")

# read analysis (x_a)
x_a = np.zeros((4,401,40))
mn  = 'im' # Choose 'pm' for perfect model, 'im' fpr imperfect model.
x_a[0,:,:] += np.genfromtxt('txt_file_new/exp1_1_x_a_'+mn+'.txt')
x_a[1,:,:] += np.genfromtxt('txt_file_new/exp1_2_x_a_'+mn+'.txt')
x_a[2,:,:] += np.genfromtxt('txt_file_new/exp1_3_x_a_'+mn+'.txt')
x_a[3,:,:] += np.genfromtxt('txt_file_new/exp1_4_x_a_'+mn+'.txt')

#%%  Calculate Error

rmse = np.zeros((401,4))
bias = np.zeros((401,4))

def RMSE(diff):
    return np.sqrt(np.sum((diff)**2)/40)

# RMSE
j = 0
tt= 0
for j in range(4):
    for tt in range(401):
        diff = x_a[j,tt,:] - x_t[tt,:]
        rmse[tt,j] = RMSE(diff)

# Mean Bias
i = 0
for i in range(4):
    diff = x_a[i,:,:] - x_t
    bias[:,i] = np.mean(diff, axis=1)

#%% Plot
'''
cn = ['g','r','b','k'] #color name
time= np.arange(401)

plt.figure(figsize=(12,7))
plt.grid()
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.8)
plt.axvline(x=50, color='r', linestyle='-', linewidth=1.8)
for i in range(4):
    plt.plot(time, rmse[:,i], color=cn[i], linestyle='-', label='$RMSE_{'+str(i+1)+'}$', marker='o', markersize=3.5)
    plt.plot(time, bias[:,i], color=cn[i], linestyle='--', label='$Bias_{'+str(i+1)+'}$')
plt.xlim([0,400])
plt.ylim([-2,7])
plt.yticks(np.arange(-1,8),fontsize=13)
plt.xlabel('t',fontsize=15)
plt.xticks(np.arange(0,401,20),0.05*np.arange(0,401,20),fontsize=13)
plt.legend(loc='upper right', numpoints=1, prop={'size':10})
plt.title('RMSE & Bias of Experiment 1',fontsize=17,fontweight='bold',loc='left')
plt.title('${Perfect\ Model}$',fontsize=16,loc='right')
#plt.title('${Imperfect\ Model}$',fontsize=16,loc='right')
plt.savefig('plot_new/verify_exp1_'+mn+'.png',dpi=300)
plt.tight_layout()
plt.show()
'''
#%% Plot zoom-in region

cn = ['g','r','b','k'] #color name
time= np.arange(81)

plt.figure(figsize=(10,7))
plt.grid()
plt.axhline(y=0, color='k', linestyle='-', linewidth=0.8)
plt.axvline(x=50, color='r', linestyle='-', linewidth=1.8)
for i in range(4):
    plt.plot(time, rmse[0:81,i], color=cn[i], linestyle='-', label='$RMSE_{'+str(i+1)+'}$', marker='o', markersize=3.5)
    plt.plot(time, bias[0:81,i], color=cn[i], linestyle='--',label='$Bias_{'+str(i+1)+'}$')
plt.xlim([0,80])
plt.ylim([-2,7])
plt.yticks(np.arange(-1,8),fontsize=13)
plt.xlabel('t',fontsize=15)
plt.xticks(np.arange(0,81,20),0.05*np.arange(0,81,20),fontsize=13)
plt.legend(loc='upper right', numpoints=1, prop={'size':10})
plt.title('RMSE & Bias of Experiment 1',fontsize=17,fontweight='bold',loc='left')
#plt.title('${Perfect\ Model}$',fontsize=16,loc='right')
plt.title('${Imperfect\ Model}$',fontsize=16,loc='right')
plt.savefig('plot_new/verify_exp1_'+mn+'_zomm_in.png',dpi=300)
plt.tight_layout()
plt.show()


