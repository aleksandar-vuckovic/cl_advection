#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 15:51:03 2019

@author: CSI\av20keco
"""

import numpy as np

k0 = -2.5
t = np.linspace(0, 1, 1000)
kappa = k0/( np.cosh(t/2)**(3/2) ) 
res = np.transpose(np.array((t, kappa)))
np.savetxt("kappa_analytic_reference.txt", res, delimiter =",")