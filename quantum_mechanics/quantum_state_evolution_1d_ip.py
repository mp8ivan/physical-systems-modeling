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
psi_inicial = np.exp(-((x-x0)**2)/(2*sigma**2))*np.exp(1j*k0*x)
norma = np.sqrt(np.sum(np.abs(psi_inicial)**2)*dx)
psi_inicial = psi_inicial/norma #la funcion de onda debe estar normalizada

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

#Representacion de la densidad para distintos instantes de tiempo:
tiempos_mostrar = [0, 0.1, 0.2, 0.5, 1.0] 
plt.figure(figsize=(10, 6))
for t_val in tiempos_mostrar:
    idx = int(round(t_val/dt))
    if idx < len(psi):
        densidad_t = np.abs(psi[idx])**2
        plt.plot(x,densidad_t,label=f't = {t_val} s')
plt.title('Evolución de la densidad de probabilidad en el pozo infinito')
plt.xlabel('x')
plt.ylabel(r'$|\psi(x,t)|^2$')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()

#Animacion:
psi_matrix = np.array(psi)
densidad = np.abs(psi_matrix)**2
max_densidad = np.max(densidad)

fig,ax = plt.subplots(figsize=(10,6))
ax.set_xlim(x[0]-0.5,x[-1]+0.5)
ax.set_ylim(0,np.max(densidad)*1.2)
ax.set_xlabel('x')
ax.set_ylabel('$|\psi (x,t)|^2$')

V_pintar = V.copy()
V_pintar[V > 1e8] = np.nan 
max_V_visible = np.nanmax(V_pintar)

if max_V_visible > 0:
    scale = (np.max(densidad)*0.8)/max_V_visible
    label_v = 'Potencial V(x) (escalado)'
else:
    scale = 1.0
    label_v = 'Potencial V(x)'

ax.plot(x,V_pintar*scale,color='gray',linestyle='--',linewidth=1,label=label_v)

if max_V_visible > 0:
    ax.fill_between(x,V_pintar*scale,np.max(densidad)*2,color='gray',alpha=0.1)
else:
    ax.axvline(x[0],color='k',linestyle='--',linewidth=1)
    ax.axvline(x[-1],color='k',linestyle='--',linewidth=1)

ax.legend(loc='upper right')
line, = ax.plot([],[],'r-',lw=2)
time_text = ax.text(0.02,0.95,'',transform=ax.transAxes)

ax.grid(True,linestyle='--',alpha=0.5)

def init():
    line.set_data([],[])
    return line,

def animate(i):
    skip = 50 
    idx = i*skip
    if idx < len(densidad):
        line.set_data(x,densidad[idx])
        time_text.set_text(f't = {t[idx]:.2f}')
    return line,time_text

ani = FuncAnimation(fig,animate,init_func=init,frames=len(t)//50,interval=20,blit=True)
plt.close()
display(HTML(ani.to_jshtml()))