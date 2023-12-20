import numpy as np
import matplotlib.pyplot as plt

L = 5.5
N = np.arange(1, 11)*1e6
print(N)

gs = [276.9, 281.7, 284.6, 268.9, 273.7, 275.1, 519.2, 505.7]
ls = [728, 19178, 25301, 46558, 81023, 140464, 188679, 161815]
r = 6e9


# size of figure
plt.figure(figsize=(10, 7))
for i in range(0, 8):
    if i==1:
        continue
    p = 2**i
    T = ((L+1)*N/p + L*N/p * (p-1)/p * gs[i] + 3*ls[i])/r
    print("\nproceccor",  f'{p}: T = ', T)
    plt.plot(N, T, '-o', label=f'{p} processors')




ax = plt.gca()
plt.ticklabel_format(axis='x',style='sci',scilimits=(6,6), useOffset=(False))
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height*0.8])

# Put a legend to the right of the current axis
plt.legend(bbox_to_anchor=(1.28, 0.82), loc = 'right')
plt.xlabel("$N$", fontsize=14)
plt.xticks(N)
plt.ylabel('expected time for SpMV (sec)', fontsize=14)
plt.show()
