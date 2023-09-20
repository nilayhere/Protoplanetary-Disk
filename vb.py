import matplotlib.pyplot as plt
import numpy as np
K_b=1.38066e-16
T=280
pi=3.141592654
rhom=3.6

from matplotlib import ticker
x=np.logspace(-4,4,500)
y=np.logspace(-4,4,500)
[X,Y]=np.meshgrid(x,y)
v_B =np.sqrt((8*K_b*T)*(rhom*(4/3)*pi*X**3+rhom*(4/3)*pi*Y**3)/(pi*(rhom*(4/3)*pi*X**3)*(rhom*(4/3)*pi*Y**3)))
plt.contourf(X,Y,v_B,locator = ticker.LogLocator() )
cbar=plt.colorbar(label="Relative velocity (cm/s)")
plt.xscale("log")
plt.title('Relative Brownian Velocity')
plt.xlabel("Size of Particle (cm)")
plt.ylabel("Size of Particle (cm)")
plt.yscale("log")
plt.show()