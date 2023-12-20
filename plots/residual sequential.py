# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 00:20:20 2023

@author: Martine
"""

import matplotlib.pyplot as plt
import numpy as np


# total time
x = np.arange(0, 127)
r  = [999.999847, 973.301432, 846.849390, 721.278538, 611.456663, 517.782291, 438.346041, 371.076274, 314.124690, 265.912840, 225.100118, 190.551308, 161.305128, 136.547743, 115.590168, 97.849195, 82.831135, 70.118073, 59.356233, 50.246138, 42.534276, 36.006044, 30.479775, 25.801688, 21.841602, 18.489316, 15.651544, 13.249319, 11.215791, 9.494373, 8.037160, 6.803603, 5.759375, 4.875416, 4.127129, 3.493690, 2.957473, 2.503555, 2.119305, 1.794031, 1.518680, 1.285590, 1.088276, 0.921245, 0.779851, 0.660158, 0.558836, 0.473065, 0.400458, 0.338995, 0.286966, 0.242922, 0.205638, 0.174076, 0.147359, 0.124742, 0.105596, 0.089389, 0.075669, 0.064056, 0.054224, 0.045902, 0.038857, 0.032893, 0.027844, 0.023571, 0.019953, 0.016891, 0.014298, 0.012104, 0.010246, 0.008673, 0.007342, 0.006215, 0.005261, 0.004454, 0.003770, 0.003192, 0.002702, 0.002287, 0.001936, 0.001639, 0.001387, 0.001174, 0.000994, 0.000842, 0.000712, 0.000603, 0.000511, 0.000432, 0.000366, 0.000310, 0.000262, 0.000222, 0.000188, 0.000159, 0.000135, 0.000114, 0.000096, 0.000082, 0.000069, 0.000059, 0.000050, 0.000042, 0.000035, 0.000030, 0.000025, 0.000022, 0.000018, 0.000015, 0.000013, 0.000011, 0.000009, 0.000008, 0.000007, 0.000006, 0.000005, 0.000004, 0.000003, 0.000003, 0.000002, 0.000002, 0.000002, 0.000001, 0.000001, 0.000001, 0.000001]

# size of figure
plt.figure(figsize=(10, 7))

# total time plot
plt.plot(x, r)


plt.xlabel("Iteration", fontsize = 15)
plt.ylabel("Residual", fontsize = 15)
# plt.xticks(r)
# plt.yscale("log")
plt.grid()
plt.ticklabel_format(axis='x',style='sci', useOffset=(False))
plt.show() 