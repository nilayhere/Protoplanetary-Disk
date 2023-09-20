from constants import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

#input parameters
R=1.49e13   # 1AU
T=280       #Temperature
sigma0=1700  # Surface density
X_m=2*10**-15 
rhom=3.6     #Solid density of particles
c_g=34300    #local thermal speed of gas
ya=1.6       #constant calculated by ormel & Cuzzi
vbt=np.zeros(((4,500,500)))
for k in range(4):
    alpha=0.00001*10**k
    #calculations
    sigma=sigma0*(R/AU)**(-3/2)
    ohm=np.sqrt(G*MSUN/R**3)
    c_s=np.sqrt(K_b*T/mump)
    H=c_s/ohm
    rhog=sigma/(H*np.sqrt(2*pi))
    v_g=np.sqrt(alpha*8/pi)*c_s
    lamda=mump/(rhog*X_m)
    v_T=alpha*c_s*H
    v_m=0.5*lamda*c_s*np.sqrt(8/pi)
    Re=v_T/v_m
    t_L=1/ohm
    t_ita=t_L/np.sqrt(Re)
    x=np.logspace(-4,4,500)
    y=np.logspace(-4,4,500)
    v_t=np.zeros((500,500))
    [X,Y]=np.meshgrid(x,y)
    t1=(rhom*X)/(rhog*c_g)
    t2=(rhom*Y)/(rhog*c_g)
    St1=t1/t_L
    St2=t2/t_L
    ep=St2/St1
    ep1=St1/St2
    for i in range(500):
        for j in range(500):
            if t1[i,j] < t_ita and t1[i,j] > t2[i,j]:
                v_t[i,j]=np.sqrt((v_g**2)*((St1[i,j]-St2[i,j])/(St1[i,j]+St2[i,j]))*((St1[i,j]**2)/(St1[i,j]+Re**-0.5)-(St2[i,j]**2)/(St2[i,j]+Re**-0.5)))
            elif t2[i,j] < t_ita and t2[i,j] > t1[i,j]:
                v_t[i,j]=np.sqrt((v_g**2)*((St2[i,j]-St1[i,j])/(St1[i,j]+St2[i,j]))*((St2[i,j]**2)/(St2[i,j]+Re**-0.5)-(St1[i,j]**2)/(St1[i,j]+Re**-0.5)))
            elif t_ita <= t1[i,j] <= t_L and t1[i,j] > t2[i,j]:
                v_t[i,j]=np.sqrt((v_g**2)*(St1[i,j])*(2*ya-1-ep[i,j]+(2/(1+ep[i,j]))*(1/(1+ya)+(ep[i,j]**3)/(ya+ep[i,j]))))
            elif t_ita <= t2[i,j] <= t_L and t2[i,j] > t1[i,j]:
                v_t[i,j]=np.sqrt((v_g**2)*(St2[i,j])*(2*ya-1-ep1[i,j]+(2/(1+ep1[i,j]))*(1/(1+ya)+(ep1[i,j]**3)/(ya+ep1[i,j]))))
            elif t1[i,j] > t_L and t1[i,j] > t2[i,j] :
                v_t[i,j]=np.sqrt((v_g**2)*(1/(1+St1[i,j])+1/(1+St2[i,j])))
            elif t2[i,j] > t_L and t2[i,j] > t1[i,j] :
                v_t[i,j]=np.sqrt((v_g**2)*(1/(1+St2[i,j])+1/(1+St1[i,j])))
vbt[k]=np.sqrt(v_t**2+(8*K_b*T)*(rhom*(4/3)*pi*X**3+rhom*(4/3)*pi*Y**3)/(pi*(rhom*(4/3)*pi*X**3)*(rhom*(4/3)*pi*Y**3)))
fig, axs = plt.subplots(2, 2)
axs[0, 0].contourf(X,Y,vbt[0],10, cmap=plt.cm.plasma)
axs[0, 0].set_title('alpha=0.00001')
plt.yscale("log")
plt.yscale("log")
axs[0, 1].contourf(X,Y,vbt[1],10, cmap=plt.cm.plasma)
axs[0, 1].set_title('alpha=0.0001')
axs[1, 0].contourf(X,Y,vbt[2],10, cmap=plt.cm.plasma)
axs[1, 0].set_title('alpha=0.001')
axs[1, 1].contourf(X,Y,vbt[3],10, cmap=plt.cm.plasma)
axs[1, 1].set_title('allpha=0.01')

for ax in axs.flat:
    ax.set(xlabel='Sizw', ylabel='y-label')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
plt.title('Relative Velocity')

plt.show()
