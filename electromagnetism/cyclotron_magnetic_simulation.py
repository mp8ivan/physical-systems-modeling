import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.mlab as mlab
from scipy.integrate import odeint
from scipy import signal
from mpl_toolkits.axes_grid1 import host_subplot
from scipy.signal import find_peaks
import random

DNI=74534379
random.seed(DNI)

dV=random.randrange(400,600,1)*10 #dV entre las dos placas desde 4000 V hasta 6000 V

print('Voltaje entre las placas: {:} V'.format(dV))

#Partículas que vamos a tratar con sus respectivas masas y cargas:
iones={'1p+': [1.67E-27, 1.6E-19],
       '1H+': [1.67E-27, 1.6E-19], '2H+': [3.32E-27, 1.6E-19],'3H+': [4.98E-27, 1.6E-19],
       '4He++': [6.688E-27, 3.2E-19],}
ni = len(iones)

def ciclotron(z,t,par): #función del ciclomotor
    
    x,vx,y,vy=z
    yy=np.abs(y)
    if yy < 0.5*h :
      cB=0
      cE=1
    else:
      cB=1
      cE=0
    r=np.sqrt((yy-0.5*h)**2+x**2) #coordenada radial
    if r > Rc:
        cE=0
        cB=0
    coswt=math.cos(omega*t) 
    s=math.copysign(1,vy) 
    cuad=signal.square(omega*t) 
    dzdt=[vx,qm*vy*cB*B,vy,qm*cE*E*s-qm*vx*cB*B]
    return dzdt

ion='4He++'
m=iones[ion][0] #masa
q=iones[ion][1] #carga

qm=q/m 
B=1  
h=0.025
Rc=0.1

E=dV/h 
omega=1.5*qm*B #frecuencia angular
print('Frecuencia angular: {:.3E} rad/s'.format(omega))

#Código para cuando el ion se dirige hacia una de las placas
tI = h*np.sqrt(2/dV/qm)
print('Tiempo en recorrer el condensador: {:.3E} s'.format(tI))
a=q*dV/m/h
print('Aceleración a la salida del condensador: {:.3E} m/s'.format(a))
v=np.sqrt(2*h*a)
print('Velocidad a la salida del condensador: {:.3E} m/s'.format(v))
R = np.sqrt(2*dV*m/(q*B**2))
print('Radio de la trayectoria curvilínea: {:.3E} m'.format(R))

#Código para cuando el ion llega a la placa
tD=np.pi/qm/B
print('T/2 = {:.3E} s'.format(tD))
print('f = {:.3E} Hz'.format(1/tD))

AE=dV*q
print('Energía que gana en cada paso AEc= {:.3E} J'.format(AE))

tf=100*tD #tiempo de simulación
par=[qm,E,B,h,Rc,omega]

nt=100000
z0=[0.0,0.0,0.0,0.0]    
t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6

z=odeint(ciclotron,z0,t,args=(par,),atol=abserr, rtol=relerr)

#Código para la gráfica del movimiento del ión:
f0 = plt.figure()
f0.set_dpi(100)
f0.set_size_inches(6,11)
f0.suptitle("Ion 1p+", fontsize=12)
ax01 = plt.subplot2grid((2, 2), (0, 0),colspan=2)
ax02 = plt.subplot2grid((2, 2), (1, 0),colspan=2)

ax01.set_xlim(-1.5*Rc,1.5*Rc)
ax01.set_ylim(-1.2*Rc,1.2*Rc)

ax02.set_xlim(0,tf)

ax02.set_xlabel("t (s)")
ax02.set_ylabel("Ec (J)")

line, = ax01.plot(z[:,0],z[:,2], linewidth=2)

Ecn=m*(z[:,1]**2+z[:,3]**2)/2
line2, = ax02.plot(t,Ecn,linewidth=2)

ax02.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax02.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))

circle1=patches.Circle((-0.7*Rc,0.3*Rc),0.05*Rc,fill=False)
circle2=patches.Circle((-0.7*Rc,0.3*Rc),0.01*Rc)
circle3=patches.Circle((-0.7*Rc,-0.3*Rc),0.05*Rc,fill=False)
circle4=patches.Circle((-0.7*Rc,-0.3*Rc),0.01*Rc)

pac1 = patches.Arc([0,0.5*h], 2*Rc, 2*Rc, angle=0, theta1=0, theta2=180)
pac2 = patches.Arc([0,-0.5*h], 2*Rc, 2*Rc, angle=180, theta1=0, theta2=180)
line1=patches.Arrow(-Rc,0.5*h,2.15*Rc,0.,width=0)
line2=patches.Arrow(-Rc,-0.5*h,2.15*Rc,0.,width=0)
circle5=patches.Circle((1.28*Rc,0),0.14*Rc,fill=False)
line3=patches.Arrow(1.17*Rc,0.03*Rc,0.05*Rc,0,width=0)
line4=patches.Arrow(1.22*Rc,0.03*Rc,0,-0.08*Rc,width=0)
line5=patches.Arrow(1.22*Rc,-0.05*Rc,0.05*Rc,0,width=0)
line6=patches.Arrow(1.27*Rc,-0.05*Rc,0,0.08*Rc,width=0)
line7=patches.Arrow(1.27*Rc,0.03*Rc,0.05*Rc,0,width=0)
line8=patches.Arrow(1.32*Rc,0.03*Rc,0,-0.08*Rc,width=0)
line9=patches.Arrow(1.32*Rc,-0.05*Rc,0.05*Rc,0,width=0)
line10=patches.Arrow(1.37*Rc,-0.05*Rc,0,0.08*Rc,width=0)
         
ax01.add_patch(circle1)
ax01.add_patch(circle2)
ax01.add_patch(circle3)
ax01.add_patch(circle4)
ax01.add_patch(circle5)

ax01.text(-0.65*Rc,0.3*Rc, "B",fontsize=20)
ax01.text(-0.65*Rc,-0.3*Rc, "B",fontsize=20)

ax01.add_patch(pac1)
ax01.add_patch(pac2)
ax01.add_patch(line1)
ax01.add_patch(line2)
ax01.add_patch(line3)
ax01.add_patch(line4)
ax01.add_patch(line5)
ax01.add_patch(line6)
ax01.add_patch(line7)
ax01.add_patch(line8)
ax01.add_patch(line9)
ax01.add_patch(line10)

rect1 = patches.Rectangle((-0.252,-0.0125),0.503,0,linewidth=1,edgecolor='black',facecolor='none')
rect2 = patches.Rectangle((-0.252,0.0146),0.503,0,linewidth=1,edgecolor='black',facecolor='none')

rect3 = patches.Rectangle((0.25,-0.0125),0.06,0,angle=-20,linewidth=0.7,edgecolor='black',facecolor='none')
rect4 = patches.Rectangle((0.25,0.014),0.06,0,angle=19,linewidth=0.7,edgecolor='black',facecolor='none')
ax01.add_patch(rect1)
ax01.add_patch(rect2)
ax01.add_patch(rect3)
ax01.add_patch(rect4)

xt = np.linspace(Rc*1.2,Rc*1.3,100)
yt = 0.02*np.sin(xt*150+2)
ax01.plot(xt, yt, 'black',linewidth= 2)

ax01.set_xlabel("x")
ax01.set_ylabel("y")

plt.show()