import matplotlib.pyplot as plt
import numpy as np
X = np.genfromtxt('xm.dat')
Y = np.genfromtxt('ym.dat')
vbt = np.genfromtxt('vbt.dat')
plt.contourf(X,Y,(vbt))
plt.colorbar()
#cbar=plt.colorbar(label="Relative velocity (cm/s)")
plt.xscale("log")
plt.title('Relative Velocity')
plt.xlabel("Size of Particle (cm)")
plt.ylabel("Size of Particle (cm)")
plt.yscale("log")
plt.show()