# -*- coding: utf-8 -*-
"""
Create full obs for NMC method
"""

import numpy as np
import matplotlib.pyplot as plt
import os

path = 'D:\\Desktop file\\All Files\\NTU課程資料\\大氣系課程\\資料同化-連國淵老師\\Project\\code'
os.chdir(path)

# read nature run (x_t) and obs error
x_t = np.genfromtxt("txt_file_new/x_t.txt")
e1  = np.genfromtxt("txt_file_new/e1.txt")  # sigma = 0.1
e2  = np.genfromtxt("txt_file_new/e2.txt")  # sigma = 0.05

full_obs_1 = x_t[0:400,:] + e1
full_obs_2 = x_t[0:400,:] + e2

np.savetxt('txt_file_new/NMC_full_obs_e1.txt', full_obs_1)
np.savetxt('txt_file_new/NMC_full_obs_e2.txt', full_obs_2)