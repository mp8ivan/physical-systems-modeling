import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from matplotlib.animation import FuncAnimation
from IPython.display import HTML,display

#Parametros:
N = 300 #numero de puntos de la simulacion
x0 = 2.5 #posicion del centro del pozo
sigma = 0.5 #anchura del pozo
k0 = 20.0 #momento inicial
t_max = 1.0 #tiempo total de simulacion

#Posicion:
x = np.linspace(0,10,N)
dx = x[1]-x[0] 

#Tiempo:
dt = 0.0001 
t = np.arange(0,t_max,dt)

#Potencial:
V = np.zeros(N)
V[0] = 1e10 
V[-1] = 1e10

#Condiciones iniciales:
En = (1**2*np.pi**2)/(10**2)
psi_inicial = np.sqrt(2/10)*np.sin(1*np.pi*x/10) #estado estacionario del pozo infinito
t_final = t[-1]
psi_analitica_final = np.sqrt(2/10)*np.sin(1*np.pi*x/10)*np.exp(-1j*En*t_final)
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

#Estado estacionario:
error_cuadratico = np.abs(psi[-1]-psi_analitica_final)**2
mse = np.mean(error_cuadratico)
print(f"Error cuadratico medio (MSE) en t={t_final}: {mse:.2e}")

#Grafica:
plt.figure(figsize=(10,5))
plt.plot(x,np.real(psi[-1]),'b-',label='Solucion numerica (parte real)')
plt.plot(x,np.real(psi_analitica_final),'r--',label='Solucion analitica (parte real)')
plt.title(f'Comparacion en t = {t_final} s (estado n=1)')
plt.xlabel('x')
plt.ylabel(r'$Re(\psi)$')
plt.legend()
plt.grid(True,alpha=0.3)
plt.show()