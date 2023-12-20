import matplotlib.pyplot as plt
import numpy as np


# total time
seq  = [3.008200,  8.667018, 14.114228, 19.233556, 24.981274, 30.560386, 36.753406, 43.135601, 49.503151, 55.327048]
p1   = [2.755350,  8.065305, 13.256462, 18.646133, 24.410600, 30.267751, 36.578800, 42.690410, 49.126685, 55.106130]
p4   = [4.826389, 10.333270, 15.927663, 21.616442, 27.243772, 33.123837, 39.195125, 45.013995, 50.901300, 56.877093]
p8   = [4.264957,  9.156307, 13.566185, 18.323920, 23.069802, 27.886642, 32.842498, 37.534052, 43.237264, 47.778738]
p16  = [2.950643,  6.042837,  9.402742, 12.615760, 15.880515, 19.057681, 22.020167, 25.480894, 29.289792, 32.325093]
p32  = [2.055643,  4.352681,  6.457889,  8.746365, 11.034820, 13.201941, 15.621388, 17.829208, 20.054626, 22.320199]
p64  = [2.012003,  4.172485,  6.362240,  8.392697, 10.630114, 12.686794, 14.968619, 17.105230, 19.107501, 21.197671]
p128 = [1.814687,  3.645621,  5.442193,  7.318216,  8.946156, 10.912390, 13.489076, 14.154535, 17.115937, 19.092282]

# plots
p = [p1, p4, p8, p16, p32, p64, p128]
N = [1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 10e6]

# size of figure
plt.figure(figsize=(10, 7))

# total time plot
for i in range(len(p)):
    if (i == 0):
        text = str(2**(i)) + " processors"
    else:
        text = str(2**(i + 1)) + " processors"
        
    plt.plot(N, p[i], '-o', label = text)

plt.plot(N, seq, 'k--', label = "sequential")

plt.xlabel("N", fontsize = 15)
plt.ylabel("total time of parallel program (sec)", fontsize = 15)
plt.legend(bbox_to_anchor=(1.22, 0.84), loc = 'right')
plt.xticks(N)
plt.ticklabel_format(axis='x',style='sci',scilimits=(6,6), useOffset=(False))
plt.show() 
    
# finding outlinks
p1   = [0.032848, 0.068505, 0.098525, 0.134754, 0.171613, 0.202818, 0.236698, 0.268677, 0.307274, 0.338632]
p4   = [0.048333, 0.090329, 0.138182, 0.179477, 0.211027, 0.270906, 0.304337, 0.342961, 0.375202, 0.417589]
p8   = [0.062522, 0.118649, 0.165981, 0.206418, 0.247830, 0.303742, 0.314455, 0.384707, 0.415487, 0.458770]
p16  = [0.067052, 0.139290, 0.175778, 0.245318, 0.296451, 0.318761, 0.467304, 0.482527, 0.497616, 0.531300]
p32  = [0.084790, 0.176081, 0.188129, 0.295942, 0.361219, 0.364205, 0.572083, 0.564202, 0.621356, 0.617086]
p64  = [0.139189, 0.234795, 0.305742, 0.350656, 0.378537, 0.577608, 0.694788, 0.703750, 0.746674, 0.768835]
p128 = [0.220669, 0.372517, 0.498510, 0.541710, 0.572907, 0.782781, 0.819860, 0.845858, 0.856934, 0.884273]

# plots
p = [p1, p4, p8, p16, p32, p64, p128]
N = [1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 10e6]

# size of figure
plt.figure(figsize=(10, 7))

# total time plot
for i in range(len(p)):
    if (i == 0):
        text = str(2**(i)) + " processors"
    else:
        text = str(2**(i + 1)) + " processors"
    plot2 = plt.plot(N, p[i], '-o', label = text)
    
plt.xlabel("N", fontsize=15)
plt.ylabel("time for finding outlinks (sec)", fontsize = 15)  
plt.legend(bbox_to_anchor=(1.22, 0.86), loc = 'right')
plt.xticks(N)
plt.ticklabel_format(axis='x',style='sci',scilimits=(6,6), useOffset=(False))
plt.show()





# time per iteration
p1   = [0.020471, 0.060111, 0.098422, 0.137610, 0.179108, 0.222311, 0.266942, 0.311711, 0.356283, 0.399783]
p4   = [0.037219, 0.078689, 0.120483, 0.162375, 0.203230, 0.247112, 0.290415, 0.333489, 0.374467, 0.418382]
p8   = [0.032733, 0.069535, 0.102381, 0.137414, 0.171857, 0.207765, 0.243244, 0.277842, 0.318001, 0.351273]
p16  = [0.022259, 0.045884, 0.070209, 0.093564, 0.117125, 0.140848, 0.160946, 0.186796, 0.213526, 0.235915]
p32  = [0.014856, 0.031524, 0.047154, 0.063415, 0.079577, 0.095836, 0.111746, 0.128336, 0.143505, 0.160197]
p64  = [0.013974, 0.029450, 0.044988, 0.060253, 0.076123, 0.090234, 0.105754, 0.121537, 0.135146, 0.150484]
p128 = [0.011353, 0.023887, 0.036586, 0.049958, 0.061643, 0.074948, 0.093059, 0.097998, 0.118965, 0.133422]

# plots
p = [p1, p4, p8, p16, p32, p64, p128]
N = [1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 10e6]

# size of figure
plt.figure(figsize=(10, 7))

# total time plot
for i in range(len(p)):
    if (i == 0):
        text = str(2**(i)) + " processors"
    else:
        text = str(2**(i + 1)) + " processors"
    plot3 = plt.plot(N, p[i], '-o', label = text)
    
plt.xlabel("N", fontsize = 15)
plt.ylabel("time per iteration (sec)", fontsize = 15)
plt.legend(bbox_to_anchor=(1.22, 0.86), loc = 'right')
plt.xticks(N)
plt.ticklabel_format(axis='x',style='sci',scilimits=(6,6), useOffset=(False))
plt.show()



