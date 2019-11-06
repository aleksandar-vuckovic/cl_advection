#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 09:24:42 2019

@author: CSI\av20keco
"""

import numpy as np
import matplotlib.pyplot as plt

lin = np.linspace(-4, 6, 100)

# Curve steepness
k = 2
L = 2
x0 = 4

logistic =  (1 / (1 + np.exp( - k * -lin )) + 1) *  (1 / (1 + np.exp( -k * (lin - x0))) + 1) -1.5

plt.plot(lin, logistic)
plt.axvline(0); plt.axvline(x0)