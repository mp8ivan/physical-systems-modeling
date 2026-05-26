import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.mlab as mlab
from scipy.integrate import odeint
from scipy import signal
import matplotlib.animation as animation
from matplotlib.pylab import *
from matplotlib import cm
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal


fig=plt.figure()
fig.set_dpi(100)

tf=15.0 #tiempo de simulacion
nt=500  #numero de intervalos de tiempo


#  Momentos de inercia

I1=4.0
I2=2.0
I3=1.0

# Valores iniciales
z1_0= 0.100                                       # WX0
z2_0= 2.000                                         # wy0 
z3_0= 0.100                                            # Wz0



par=[I1,I2,I3]
print('Momentos de inercia ',par)




omel= np.sqrt((I1-I2)*(I1-I3)/(I2*I3))*z1_0 #frecuencia angular de precesion del nucleo
timel= 2*np.pi/omel #periodo de precesion del nucleo
print(' precesion nucleo=', omel, timel)




lo2= (I1*z1_0)**2 + (I2*z2_0)**2 + (I3*z3_0)**2
eoci= 0.5*(  I1*z1_0**2 + I2*z2_0**2 + I3*z3_0**2  )





# Definiendo funcion

def rotacion(z,t,par):
    z1,z2,z3 = z
    dzdt = [(I2-I3)*z2*z3/I1,(I3-I1)*z1*z3/I2,(I1-I2)*z1*z2/I3]
    return dzdt



# Llamada a odeint que resuelve las ecuaciones de movimiento


z0=[z1_0,z2_0,z3_0] #Valores iniciales   

print('Estado inicial',z0)

t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6

z=odeint(rotacion,z0,t,args=(par,),atol=abserr, rtol=relerr)

plt.close('all')

# Extrae las velocidades angulares
velocidad_angular_x = z[:, 0]
velocidad_angular_y = z[:, 1]
velocidad_angular_z = z[:, 2]

# Usa scipy.signal.find_peaks para encontrar los picos
picosX, _ = signal.find_peaks(velocidad_angular_x)

plt.figure()
plt.xlabel('t (s)')
plt.ylabel('$\omega_x$ (rad/s)')
plt.plot(t, velocidad_angular_x,linewidth=1,color='blue')
plt.plot(t[picosX], velocidad_angular_x[picosX],'*',color='r')


picosY, _ = signal.find_peaks(velocidad_angular_y)

plt.figure()
plt.xlabel('t (s)')
plt.ylabel('$\omega_y$ (rad/s)')
plt.plot(t, velocidad_angular_y,linewidth=1,color='blue')
plt.plot(t[picosY], velocidad_angular_y[picosY],'*',color='r')


picosZ, _ = signal.find_peaks(velocidad_angular_z)

plt.figure()
plt.xlabel('t (s)')
plt.ylabel('$\omega_z$ (rad/s)')
plt.plot(t, velocidad_angular_z,linewidth=1,color='blue')
plt.plot(t[picosZ], velocidad_angular_z[picosZ],'*',color='r')



plt.figure()
plt.xlabel('t (s)')
plt.ylabel('$\omega$ (rad/s)')
plt.plot(t,velocidad_angular_x,label='$\omega_x$(t)',linewidth=1,color='orange')
plt.plot(t,velocidad_angular_y,label='$\omega_y$(t)',linewidth=1,color='black')
plt.plot(t,velocidad_angular_z,label='$\omega_z$(t)',linewidth=1,color='blue')
plt.legend(loc='lower right')

# Calculo momento angular y energía cinética

ex=z[:,0]
ey=z[:,1]
ez=z[:,2]

lx=I1*ex
ly=I2*ey
lz=I3*ez

l22= np.sqrt(lx**2+ly**2+lz**2)

plt.figure()
plt.xlabel('t (s)')
plt.ylabel('L (kg·$\mathrm{m^2}$)')
plt.plot(t[:],l22[:],linewidth=1 )
#plt.ylim(4.016,4.026)

enecine=0.5*(I1*z[:,0]**2+I2*z[:,1]**2+I3*z[:,2]**2)


plt.figure()
plt.xlabel('t (s)')
plt.ylabel('$\mathrm{E_c}$ (J)')
plt.plot(t[:],enecine[:],color='orange',linewidth=1)
#plt.ylim(4.020,4.030)


############## Recorrido de la velocidad angular y el momento angular


wxx=np.linspace(0,z1_0,10)
wxy=np.linspace(0,z2_0,10)
wxz=np.linspace(0,z3_0,10)


fig= plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')
ax.view_init(10,-40)
ax.plot3D(ex,ey,ez,color='orange')
ax.plot3D(lx,ly,lz,color='blue')



loxx=np.linspace(0,I1*z1_0,10)
loxy=np.linspace(0,I2*z2_0,10)
loxz=np.linspace(0,I3*z3_0,10)

fig= plt.figure()
ax = fig.add_subplot(1,1,1,projection='3d')
ax.view_init(10,-40)
ax.plot3D(lx,ly,lz)
