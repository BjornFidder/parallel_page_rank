import numpy as np
import matplotlib.pyplot as plt

L = 5.5
N = np.arange(11)*1e6

gs = [276.9, 281.7, 284.6, 268.9, 273.7, 275.1, 519.2, 505.7]
ls = [728, 19178, 25301, 46558, 81023, 140464, 188679, 161815]
r = 6e9

plt.figure()
for i in range(0, 8):
    if i==1:
        continue
    p = 2**i
    T = ((L+1)*N/p + L*N/p * (p-1)/p * gs[i] + 3*ls[i])/r
    plt.plot(N, T, label=f'{p} processors')

ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height*0.8])

# Put a legend to the right of the current axis
ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)


plt.xlabel("$N$")
plt.ylabel('expected time for SpMV (sec)')
plt.show()
