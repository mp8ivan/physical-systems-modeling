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
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema


fig=plt.figure()
fig.set_dpi(100)
fig.set_size_inches(7,6.5)


# Datos a modificar en la simulacion 

m=500.0                           # Valor de las masa
mati=5.9*10**(24)
gr=6.67*10**(-11)
radi= 6378000  # radio tierra
tf=9500    #tiempo de simulacion
f = 1 #factor modifica
k = 1

# Para dibujar el radio de la Tierra

rapla = np.zeros(360,float)
radix = np.zeros(360,float)
radiy = np.zeros(360,float)

for j in range (0,360):
  rapla[j] = 2*j*np.pi/360
  radix[j]= radi*np.cos(rapla[j])  #  
  radiy[j]= radi*np.sin(rapla[j])  # 


par=[gr,mati,radi]


# Definiendo  gravedad

def tiro(z,t,par):
    z1,z2,z3,z4=z  #z1 = posición x; z2 = posición y; z3 = velocidad x; z4 = velocidad y
    dzdt=[z3,z4, 
          -k*z1/m,
         -k*z2/m]
    return dzdt
# Llamada a odeint que resuelve las ecuaciones de movimiento

nt=75000  #numero de intervalos de tiempo
dt=tf/nt

# Valores iniciales
z1_0= radi+300000                 # x punto  inicial
z2_0=0.0                         # y punto  inicial
z3_0= 0                          # Velocidad x 
z4_0= np.sqrt(gr*mati/z1_0)*f      # Velocidad y

z0=[z1_0,z2_0,z3_0,z4_0] #Valores iniciales   
 
t=np.linspace(0,tf,nt)
at=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6

z=odeint(tiro,z0,t,args=(par,),atol=abserr, rtol=relerr)

# Calculo del momento angular inicial y energia total inicial
# SOLO SI V TIENE SOLO COMPONENTE Y LA POSICIÓN COMPONENTE X 

lo= m*z4_0*z1_0 #momento angular 
E0=0.5*m*z4_0**2-gr*mati*m/z1_0 #energía

print('Energía = {}'.format(E0))
print('Momento angular = {}'.format(lo))
print('Velocidad y = {}'.format(z4_0))
    
# definimos el tramo radial sobre el que representamos el 
# potencial y la energía. Hay que ajustarlo a mano.

rix= np.linspace(z1_0*0.9,1.6*z1_0,100)

Pote= 0.5*lo**2/m**2/rix**2-gr*mati/rix
ene1=np.linspace(E0,E0,100)

mu = (m*mati)/(m+mati)

### Gráficas

# Recordemos z[:,0] es x 
# Recordemos z[:,1] es y
# Recordemos z[:,2] es Vx
# Recordemos z[:,3] es Vy

plt.close('all')

# Trayectoria xy
plt.figure(figsize=(4,4))
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.plot(z[:,0],z[:,1],linewidth=1)
plt.plot(radix,radiy, linewidth=1)  # dibuja el contorno de la Tierra
plt.plot(0,0,'*')

# Evolución del radio con el tiempo

rad=np.sqrt( z[:,0]**2 + z[:,1]**2)
minimos_indices = argrelextrema(rad, np.less)[0]

plt.figure(figsize=(4,4))
plt.xlabel('t (s)')
plt.ylabel('r (m)')
plt.plot(t[:],rad[:],linewidth=1)
plt.scatter(t[minimos_indices], rad[minimos_indices], color='red')

# Velocidades

plt.figure(figsize=(4,4))
plt.xlabel('v (eje x)')
plt.ylabel('v (eje y)')
plt.plot(z[:,2],z[:,3],linewidth=1)

# Evolución de la velocidad con el tiempo

vel=np.sqrt( z[:,2]**2 + z[:,3]**2)

plt.figure(figsize=(4,4))
plt.xlabel('t (s)')
plt.ylabel('v (m/s)')
plt.plot(t[:],vel[:],linewidth=1)
    
# Evolucion de la energia en funcion del tiempo     

plt.figure(figsize=(4,4))
plt.xlabel('r (m)')
plt.ylabel('E (J)')
plt.plot(rix,Pote,linewidth=0.5 )
plt.plot(rix,ene1 ,linewidth=1.5)  

#Animación
x=z[:,0]
y=z[:,1]

fig, ax = plt.subplots()
ax.set_xlim(min(x)-1e7, max(x)+1e7)
ax.set_ylim(min(y)-1e7, max(y)+1e7)
ax.plot(0, 0, "*r")
line, = ax.plot(x[:1], y[:1])
plt.plot(radix,radiy, linewidth=1) 
a = 100

ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

def Orb(i):
    
    line.set_data(x[:a*i+1], y[:a*i+1])
    
    return line

ani = animation.FuncAnimation(fig, Orb, frames=len(x)//a, interval=5)


#ani.save("orbeli.gif")

# Encuentra los puntos donde la velocidad radial es cero
radial_velocity_zeros = np.where(np.diff(np.sign(z[:, 2])))[0]

# Calcula las distancias correspondientes a estos puntos
apogee_distance = max(rad[radial_velocity_zeros])
perigee_distance = min(rad[radial_velocity_zeros])

print('Apogeo: {} metros'.format(apogee_distance))
print('Perigeo: {} metros'.format(perigee_distance))