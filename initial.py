# -*- coding: utf-8 -*-
"""
Create initial condition for DA experiment
Save:
  x_a_init.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from setting_perfect_model import *

# settings of spin-up
sigma_x0 = 0.1  # size of initial perturpation
Tspinup = 100.  # initial spin-up time

# spin-up from a random initail value
x_a_0 = sigma_x0 * np.random.randn(N)

solver = ode(lorenz96.f).set_integrator('dopri5', nsteps=10000)
solver.set_initial_value(x_a_0, 0.).set_f_params(F)
solver.integrate(Tspinup)
x_a_init = np.array(solver.y, dtype='f8')

# save the initial condition for DA experiment
np.savetxt('txt_file_new/x_a_init_perfect.txt', x_a_init)
