import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from matplotlib.animation import FuncAnimation
from IPython.display import HTML,display

#Parametros:
N = 300 #numero de puntos de la simulacion
x0 = 2.0 #posicion del centro del pozo
sigma = 0.5 #anchura del pozo
k0 = 0.0 #momento inicial
t_max = 1.0 #tiempo total de simulacion

#Posicion:
x = np.linspace(-5,5,N)
dx = x[1]-x[0] 

#Tiempo:
dt = 0.0001 
t = np.arange(0,t_max,dt)

#Potencial:
V = 0.5*(5.0**2)*x**2 #omega = 5.0, m = 1.0
V[0] = 1e10 
V[-1] = 1e10

#Condiciones iniciales:
psi_inicial = np.exp(-((x-x0)**2)/(2*sigma**2))*np.exp(1j*k0*x) #paquete gaussiano
norma = np.sqrt(np.sum(np.abs(psi_inicial)**2) * dx)
psi_inicial = psi_inicial/norma #la funcion de onda debe estar normalizada
psi_inicial[0] = 0.0
psi_inicial[-1] = 0.0

#Matrices difusion:
H = diags([-1.0/(dx**2)*np.ones(N-1),2.0/(dx**2)+V,-1.0/(dx**2)*np.ones(N-1)],offsets=[-1,0,1],shape=(N,N),format='csc')
I = sparse.eye(N,N,dtype=complex,format='csc') #matriz identidad
B = I+(1j*dt/2)*H 
M = I-(1j*dt/2)*H 

#Condiciones de contorno:
B = B.tolil() #pasamos B a lista para poder añadir las condiciones de contorno facilmente
B[0,:] = 0.0
B[0,0] = 1.0       
B[-1,:] = 0.0
B[-1,-1] = 1.0
B = B.tocsc() #B vuelve a estar en formato csc para poder usar spsolve (y para mayor eficiencia)

M = M.tolil() #pasamos M a lista para poder añadir las condiciones de contorno facilmente
M[0,:] = 0.0   
M[-1,:] = 0.0
M = M.tocsc() #M vuelve a estar en formato csc para poder usar spsolve (y para mayor eficiencia)

#Solucion:
psi_actual = psi_inicial.copy() 
psi = [psi_actual.copy()] 
for i in range(1,len(t)):
    s = M@psi_actual
    psi_actual = spsolve(B,s) #B*psi^{N+1} = s
    psi.append(psi_actual.copy())

#Potencial armonico:
k_resorte = 5.0**2           
m_num = 0.5                 
omega_fisica = np.sqrt(k_resorte / m_num)
x_cm_num = [np.sum(x*np.abs(p)**2)*dx for p in psi] #trayectoria del centro de masa
x_cl_t = x0*np.cos(omega_fisica*t) #x(t) = x0*cos(omega*t)
energias = [np.real(np.conj(p)@(H@p))*dx for p in psi]

#Graficas:
plt.figure(figsize=(10,5))
plt.plot(t,x_cm_num,'b-',label='Centro de masa numerico')
plt.plot(t,x_cl_t,'r--',label='Trayectoria clasica')
plt.title('Oscilacion en potencial armonico')
plt.xlabel('Tiempo (s)')
plt.ylabel('Posicion')
plt.legend()
plt.grid(True,alpha=0.3)
plt.show()

plt.figure(figsize=(10,5))
plt.plot(t,energias,color='green')
plt.ylim(np.mean(energias)*0.9,np.mean(energias)*1.1)
plt.title(fr'Conservacion del valor esperado de la energia $\langle E \rangle$')
plt.ylabel('Energia')
plt.xlabel('Tiempo (s)')
plt.grid(True,alpha=0.3)
plt.show()