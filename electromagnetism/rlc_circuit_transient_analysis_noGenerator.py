import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.mlab as mlab
from scipy.integrate import odeint
from scipy.signal import find_peaks

DNI=74534379
random.seed(DNI)

C = random.randrange(50,101,1)*1E-6 # entre 50 y 100 μF
L = random.randrange(10,101,1)*1E-2 # entre 0.1 y 1 H
Q_inicial = random.randrange(10,101,1)*1E-7 # entre 1 y 10  μC

print('C = {:.2E} F, L = {:.2E} H y Q_inicial = {:.2E} C'.format(C,L,Q_inicial))

#Código para el condensador RLC:
def circRCL(z,t,par): #ecuaciones del sistema
    Q_inicial,I = z  
    dzdt=[I,(-R*I - Q_inicial/C)/L]
    return dzdt

R = 10.0  #resistencia (ohm)
tf = 0.1 #tiempo de simulación

par = [R,C,L]

nt = 10000
z_inicial = [Q_inicial,0.0]    
t = np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z = odeint(circRCL,z_inicial,t,args = (par,),atol = abserr,rtol = relerr)

#Código para el gráfico de la carga en función del tiempo:
fig = plt.figure()
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf))
line1, = ax.plot(t,z[:,0],'--',label = 'Q(C)',linewidth = 2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("Q(C)")
ax.set_xlabel("t(s)")
plt.title('Carga en función del tiempo')
plt.show()

#Código para el oscilador no amortiguado:
R = 0.0  #resistencia (ohm)
tf = 0.3 #tiempo de simulación

par = [R,C,L]

nt = 10000
z_inicial = [Q_inicial,0.0]    
t = np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z = odeint(circRCL,z_inicial,t,args = (par,),atol = abserr,rtol = relerr)

freq_eq = 1/2*np.pi*np.sqrt(L*C) #valor analítico de la frecuencia
print('f_analítica = '+str(freq_eq))
print('Periodo analítico = '+str(1/freq_eq))

#Picos de la función definida anteriormente:
Qpeaks, _ = find_peaks(z[:,0])
Ipeaks, _ = find_peaks(z[:,1])

Qper = (t[Qpeaks[-1]]-t[Qpeaks[0]])/(len(Qpeaks)-1)
print('f_numérica de Q(t) = '+str(1/Qper))
Iper = (t[Ipeaks[-1]]-t[Ipeaks[0]])/(len(Ipeaks)-1)
print('f_analítica de I(t) = '+str(1/Iper))

#Código para el gráfico de la carga en función del tiempo:
fig = plt.figure('Q=f(t)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf))
line1, = ax.plot(t,z[:,0],'--',label = 'Q(C)',linewidth=2)
line1p, = ax.plot(t[Qpeaks],z[Qpeaks,0],'x')
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("Q(C)")
ax.set_xlabel("t(s)")
plt.title('Q=f(t)')
plt.show()

#Código para el gráfico de la corriente en función del tiempo:
fig = plt.figure('I=f(t)')
fig.set_dpi(100)
ax=plt.axes(xlim = (0,tf))
line2, = ax.plot(t,z[:,1],'--',label = 'I(A)', linewidth = 2)
line2p, = ax.plot(t[Ipeaks],z[Ipeaks,1],'x')
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style='sci',scilimits = (-2,2))
ax.set_ylabel("I(A)")
ax.set_xlabel("t(s)")
plt.title('I=f(t)')
plt.show()

Um= (L*(z[:,1]**2))/2 #energía de la bobina
Ue= (z[:,0]**2)/(2*C) #energía del condensador
Ut = Um + Ue #energía total

#Picos de la función de la energía:
Upeaks, _ = find_peaks(Um)
Um_per = (t[Upeaks[-1]]-t[Upeaks[0]])/(len(Upeaks)-1)
print('Frecuencia de Um(t) = '+str(1/Um_per)+'Hz')

#Código para el gráfico de las energías en función del tiempo:
fig = plt.figure('Um(J),Ue(J),Ut(J)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf))
line1, = ax.plot(t[:],Um,'--',label = 'Um(J)',linewidth = 2)
line2, = ax.plot(t[:],Ue,'--',label = 'Ue(J)',linewidth = 2)
line3, = ax.plot(t[:],Ut,'-',label = 'Ut(J)',linewidth = 2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci', scilimits = (-2,2))
ax.set_ylabel("U(J)")
ax.set_xlabel("t(s)")
plt.title('Energías en función del tiempo')
plt.show() 

#Código para el gráfico de la energía total en función del tiempo:
fig = plt.figure('Ut(J)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf))
line3, = ax.plot(t[:],Ut,'-',label = 'Ut(J)',linewidth = 2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("Ut(J)")
ax.set_xlabel("t(s)")
plt.title('Energía total en función del tiempo')
plt.show()

Et = (Q_inicial**2)/(2*C) #energía total analítica
print('Energía total analítica = {:.3E} J'.format(Et))

#Código para la oscilación subamortiguada:
R=0.3*np.sqrt(4*L/C)
print('R = {:.0f} ohm'.format(R))

par = [R,C,L]

tf = 15.0E-2 #tiempo de simulación
nt = 10000
z_inicial = [Q_inicial,0.0]    
t = np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z = odeint(circRCL,z_inicial,t,args = (par,),atol = abserr,rtol = relerr)

#Código para la gráfica de la carga en función del tiempo:
fig = plt.figure('Q(C)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf))
line1, = ax.plot(t[:],z[:,0],'--',label = 'Q(s)',linewidth = 2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("Q(C)")
ax.set_xlabel("t(s)")
plt.title('Carga en función del tiempo')
plt.show()

#Código para la gráfica de la intensidad en función del tiempo:
fig = plt.figure('I(C/s)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf))
line2, = ax.plot(t[:],z[:,1],'--',label = 'I', linewidth = 2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("I(C/s)")
ax.set_xlabel("t(s)")
plt.title('Intensidad en función del tiempo')
plt.show()

Um= (L*(z[:,1]**2))/2 #energía de la bobina
Ue= (z[:,0]**2)/(2*C) #energía del condensador
Ut = Um + Ue #energía total

#Código para las energías en función del tiempo:
fig = plt.figure('Um(J),Ue(J),Ut(J)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf/2))
line1, = ax.plot(t[:],Um,'--',label='Um(J)', linewidth=2)
line2, = ax.plot(t[:],Ue,'--',label='Ue(J)', linewidth=2)
line3, = ax.plot(t[:],Ut,'-',label='Ut(J)', linewidth=2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("U(J)")
ax.set_xlabel("t(s)")
plt.show()

#Gráfica para la energía total en función del tiempo:
fig=plt.figure('Ut(J)')
fig.set_dpi(100)
ax = plt.axes(xlim = (0,tf/2))
line3, = ax.plot(t[:],Ut,'-',label = 'Ut(J)',linewidth=2)
ax.legend(loc = 'lower right')
ax.ticklabel_format(axis = 'y',style = 'sci',scilimits = (-2,2))
ax.set_ylabel("U(J)")
ax.set_xlabel("t(J)")
plt.show()

#Código para la oscilación sobreamortiguada:
R=1.5*np.sqrt(4*L/C)
print('R = {:.1f} ohm'.format(R))

par=[R,C,L]

tf=0.06 #tiempo de simulación
nt=10000
z_inicial=[Q_inicial,0.0]    
t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z=odeint(circRCL,z_inicial,t,args=(par,),atol=abserr, rtol=relerr)

#Código para el gráfico de la carga:
fig=plt.figure('Q=f(t)')
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line1, = ax.plot(t,z[:,0],'--',label='Q(C)', linewidth=2)
ax.legend(loc='lower right')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("Q(C)")
ax.set_xlabel("t(s)")
plt.show()

#Código para la gráfica de la corriente:
fig=plt.figure('I=f(t)')
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line2, = ax.plot(t[:],z[:,1],'--',label='I(A)', linewidth=2)
ax.legend(loc='lower right')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("I(A)")
ax.set_xlabel("t(s)")
plt.show()

Um= (L*(z[:,1]**2))/2 #energía de la bobina
Ue= (z[:,0]**2)/(2*C) #energía del condensador
Ut = Um + Ue #energía total

#Código para la gráfica de las energías en función del tiempo:
fig=plt.figure('Um(J),Ue(J),Ut(J)')
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line1, = ax.plot(t[:],Um,'--',label='Um(J)', linewidth=2)
line2, = ax.plot(t[:],Ue,'--',label='Ue(J)', linewidth=2)
line3, = ax.plot(t[:],Ut,'-',label='Ut(J)', linewidth=2)
ax.legend(loc='lower right')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("U(J)")
ax.set_xlabel("t(s)")
plt.show()

fig=plt.figure('Enegía total')
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf/2))
line3, = ax.plot(t,Ut,'-',label='Ut(J)', linewidth=2)
ax.legend(loc='lower right')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("Ut(J)")
ax.set_xlabel("t(s)")
plt.show()

#Código para el apartado de la energía disipada por la resistencia:
def dis_circRCL(z,t,par): #función de la energía disipada
    Q,I,Ur=z  
    dzdt=[I,-(R*I+(Q/C))/L,(I**2)*R]
    return dzdt

tf=0.1 #tiempo de simulación
par=[R,C,L]

nt=10000
z_inicial=[Q_inicial,0.0,0.0]    
t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z=odeint(dis_circRCL,z_inicial,t,args=(par,),atol=abserr, rtol=relerr)

Ue = (z[:,0]**2)/(2*C) #energía del condensador

print('Energía disipada: {:.3E} J'.format(z[:,2][-1]))

#Código para el gráfico de la energía:
fig=plt.figure('U=f(t)')
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line3, = ax.plot(t,z[:,2],'--',label='Ur(J)', linewidth=2)
line2, = ax.plot(t,Ue,'--',label='Ue(J)', linewidth=2)
ax.legend(loc='lower right')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("U(J)")
ax.set_xlabel("t(s)")
plt.show()

#Código para el gráfico de la corriente: 
fig=plt.figure('I=f(t)')
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line3, = ax.plot(t,z[:,1],'--',label='Ur(J)', linewidth=2)
ax.legend(loc='lower right')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("I(A)")
ax.set_xlabel("t(s)")
plt.show()