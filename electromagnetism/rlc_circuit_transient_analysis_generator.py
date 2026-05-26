import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.signal import find_peaks

DNI=74534379

random.seed(DNI)
R=random.randrange(100,501,1) #resistencia comprendida entre 100 ohm y 500 ohm
C=random.randrange(50,101,1)*1E-6 #capacidad comprendida entre 50 μF y 100 μF
L=random.randrange(10,51,1)*1E-2 #inductancia comprendida entre 0.1 H y 0.5  H
Ve=random.randrange(110,221,1) #voltaje eficaz comprendida entre 110 V y 220 V
f=random.randrange(50,101,1) #frecuencia comprendida entre 50 Hz y 100  Hz 

print('R = {:} ohm, C = {:.2E} F, L = {:.2E} H, Ve = {:} V y f = {:} Hz'.format(R,C,L,Ve,f))

#Cálculo analítico:
XL=2*np.pi*f*L  #reactancia inductiva
XC=1/(2*np.pi*f*C) #reactancia capacitativa
print('Reactancia inductiva: {:.1f} ohm'.format(XL))
print('Reactancia capacitativa: {:.1f} ohm'.format(XC))

Z=(R**2+(XL-XC)**2)**(1/2) #impedancia
print('Impedancia: {:.1f} ohm'.format(Z))

desfase = np.arctan((XL-XC)/R) #desfase 
print('Desfase: {:.3f} rad'.format(desfase))
print('Desfase: {:.3f} grados'.format(desfase*180/np.pi))

V0 = Ve*((2)**(1/2))
Ip= V0/Z #intensidad pico
print('Intensidad pico: {:.3} A'.format(Ip))
Ie= Ip/((2)**(1/2)) #intensidad eficaz
print('Intensidad eficaz: {:.3E} A'.format(Ie))

#Código para el oscilador amortiguado y forzado:
def circRCL(z,t,par): #función de las ecuaciones del movimiento
    Q,I=z  
    dzdt=[I,(V0*np.cos(omega*t)-Q/C-R*I)/L]
    return dzdt

V_inicial = Ve*((2)**(1/2)) #voltaje pico
omega=2.*np.pi*f  #frecuencia angular

par=[R,C,L,V0,omega]

tf=0.1 #tiempo de simulación
nt=100000 
z_inicial=[0.0,0.0]  
t=np.linspace(0,tf,nt)
abserr=1.0e-8
relerr=1.0e-6
z=odeint(circRCL,z_inicial,t,args=(par,),atol=abserr, rtol=relerr)

Vgen=V_inicial*np.cos(omega*t) #voltaje en función del tiempo

#Picos de la intensidad y el voltaje en función del teimpo:
Ipeaks,_ = find_peaks(z[:,1]) 
Vpeaks,_ = find_peaks(Vgen) 

#Cálculo numérico de la intensidad:
Ie=z[Ipeaks[-1],1]/np.sqrt(2)
print('Intensidad máxima: {:.3E} A'.format(z[Ipeaks[-1],1]))
print('Intensidad eficaz: {:.3E} A'.format(Ie))

#Cálculo numérico del periodo y la frecuencia de la intensidad:
Iper=(t[Ipeaks[-1]]-t[Ipeaks[-5]])/4 
print('Periodo de I(t) = {:.3E} s'.format(Iper))
print('Frecuencia de I(t) = {:.3E} Hz'.format(1/Iper))

#Cálculo numérico del desfase entre la intensidad y el voltaje:
desfase=2*np.pi*(t[Ipeaks[-2]]-t[Vpeaks[-2]])/Iper
print('Desfase: {:.3E} rad'.format(desfase))
print('Desfase: {:.3E} grados'.format(180/np.pi*desfase))

#Código para la gráfica del voltaje y la corriente en función del tiempo (con ejes de coordenadas distintos):
fig, ax1 = plt.subplots()
fig.set_dpi(150)
ax1.set_xlabel('t(s)')
ax1.set_ylabel('I(A)', color='b')
ax1.tick_params('y', colors='b')
ax2 = ax1.twinx()
ax1.set_xlim(xmin=0.0,xmax=tf)
line1, = ax1.plot(t,z[:,1],'-', label='I(A)', linewidth=2, color='b')
line2, = ax2.plot(t,Vgen,'--', label='V(V)', linewidth=2, color='r')
ax2.set_ylabel('V(V)', color='r')
ax2.tick_params('y', colors='r')
fig.tight_layout()
plt.show()

#Código para la gráfica de la corriente en función del tiempo:
fig=plt.figure()
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line1, = ax.plot(t,z[:,1],'--', linewidth=2)
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("I(A)")
ax.set_xlabel("t(s)")
plt.show()

#Código para la gráfica de la carga en función del tiempo:
fig=plt.figure()
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line1, = ax.plot(t,z[:,0],'--', linewidth=2)
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("Q(C)")
ax.set_xlabel("t(s)")
plt.show()

#Códgio para la gráfica del voltaje en función del tiempo:
fig=plt.figure()
fig.set_dpi(100)
ax=plt.axes(xlim=(0,tf))
line1, = ax.plot(t,Vgen,'--', linewidth=2)
ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.set_ylabel("V(V)")
ax.set_xlabel("t(s)")
plt.show()

#Código para el estudio de la resonancia:
R = 100
f_res =1/(2*np.pi*(L*C)**(1/2)) #frecuencia de resonancia
print('Frecuencia de resonancia obtenida analíticamente: {:.2f} Hz'.format(f_res))

XLf= 2*np.pi*f_res*L #reactancia inductiva
XCf= 1/(2*np.pi*f_res*C) #reactancia capacitativa
Zf = np.sqrt((R**2)+(XLf-XCf)**2) #impedancia  
Ipf = Ve*np.sqrt(2)/Zf #intensidad pico
#Valores analíticos:
print('Reactancia inductiva: {:.1f} ohm'.format(XLf))
print('Reactancia capacitativa: {:.1f} ohm'.format(XCf))
print('Impedancia: {:.1f} ohm'.format(Zf))
print('Intensidad pico: {:.1f} A'.format(Ipf))

#Cálculo de la intensidad pico en función de la frecuencia:
nf=10000
fi=np.linspace(10,100,nf)
vI=[]
for i in fi:
  XL=2*np.pi*i*L
  XC=1/(2*np.pi*i*C)
  I=Ve*np.sqrt(2)/(R**2+(XL-XC)**2)**0.5
  vI.append(I)

#Picos de la intensidad en función del tiempo:
vIa=np.array(vI) 
vIpeaks,_ = find_peaks(vIa) 
print('Intensidad pico: {:.4f} A'.format(vIa[vIpeaks][-1]))
print('Intensidad eficaz: {:.4f} A'.format(vIa[vIpeaks][-1]/np.sqrt(2)))
print('Frecuencia: {:.2f} Hz'.format(fi[vIpeaks][-1]))

#Código para la gráfica de la intensidad en función de la frecuencia:
plt.figure()
plt.plot(fi,vIa,'-b',label='Intensidad pico')
plt.plot(fi,vIa/np.sqrt(2),'-r',label='Intensidad eficaz')
plt.plot(fi[vIpeaks],vIa[vIpeaks][-1],'x')
plt.plot(fi[vIpeaks],vIa[vIpeaks]/np.sqrt(2),'o')
plt.annotate('$máx$.',xy=(fi[vIpeaks],vIa[vIpeaks]),arrowprops=dict(arrowstyle='->'),xytext=(fi[vIpeaks]*1.3,vIa[vIpeaks]))
plt.xlabel('f(Hz)')
plt.ylabel('I(A)')
plt.show()

omega=2.*np.pi*f_res #frecuencia angular
tf=0.3 #tiempo final la simulación
par=[R,C,L,V0,omega]
nt=100000 
z_inicial=[0.0,0.0] 
t=np.linspace(0,tf,nt)
abserr = 1.0e-8
relerr = 1.0e-6
z=odeint(circRCL,z_inicial,t,args=(par,),atol=abserr, rtol=relerr)
Vm=V0*np.cos(omega*t[:])

#Código para el gráfico de la intensidad y el voltaje en función del tiempo
fig=plt.figure()
fig.set_dpi(100)
fig, ax1 = plt.subplots()
ax1.set_xlabel('t(s)')
ax1.set_ylabel('I(A)', color='r')
ax1.tick_params('y', colors='r')
ax2 = ax1.twinx()
ax1.set_xlim(xmin=0,xmax=tf)
line1, = ax1.plot(t[:],z[:,1],'-', label='I(A)', linewidth=2, color='b')
line2, = ax2.plot(t[:],Vm,'--', label='V(V)', linewidth=2, color='r')
ax2.set_ylabel('V(V)', color='black')
ax2.tick_params('y', colors='black')
fig.tight_layout()
plt.show()

#Código para la gráfica de la intensidad (I=V/R) en función del tiempo:
fig=plt.figure()
fig.set_dpi(100)
fig, ax1 = plt.subplots()
ax1.set_xlabel('t(s)')
ax1.set_ylabel('I(A)', color='b')
ax1.tick_params('y', colors='b')
ax2 = ax1.twinx()
ax1.set_xlim(xmin=0,xmax=tf) #limites del eje x
line1, = ax1.plot(t[:],z[:,1], label='I(A)', linewidth=2, color='b')
line2, = ax2.plot(t[:],Vm/R, label='I=V/R', linewidth=2, color='g')
ax2.set_ylabel('V/I (A)', color='y')
ax2.tick_params('y', colors='y')
fig.tight_layout()
plt.show()