import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt(r"C:\Users\bjorn\Parallel\PageRank\plots\p128.txt", delimiter=' ').transpose()

h = data[0]
t = data[1]

g = 505.7
T0 = 175453

plt.figure()
plt.plot(h, t, '+', label='Measured data')
plt.plot(h, g*h+T0, '--', label='Least squares fit')
plt.ylim((0, 3e6))
plt.xlabel('$h$')
plt.ylabel('Time (in flop units)')
plt.legend()
plt.show()