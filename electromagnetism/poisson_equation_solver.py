import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
charge_colors = {True: '#aa0000', False: '#0000aa'}

N=21
q = 3.
L0 = -2.0
L = 4.0
h=L/(N-1)
error=0.00001 
iteraciones=0
Solucion=0
V=np.zeros((N,N))
V_temp=np.zeros((N,N))

def Pos(r):
    i=int((r[0]-L0)/h)
    j=int((r[1]-L0)/h)
    return i,j

cargas=[]
cargas.append((q,(-0.5,0)))
cargas.append((-q,(0.5,0)))

rho = np.zeros((N,N))

for carga in cargas:
	rho[Pos(carga[1])]=carga[0]

while Solucion<N**2:
    Solucion=0
    iteraciones+=1
    for i in range(1,N-1):
        for j in range(1,N-1):
            V[i,j]=1/4*(V_temp[i+1,j]+V_temp[i-1,j]+V_temp[i,j+1]+V_temp[i,j-1]+rho[i,j]*h**2)
    for i in range(N):
        for j in range(N):
            if abs(V[i,j]-V_temp[i,j])>error: 
                Solucion=0                    
            else:
                Solucion+=1 
    V_temp=np.copy(V)

print('{} Iteraciones'.format(iteraciones))

Ey , Ex = np.gradient(-np.transpose(V_temp))
x=np.linspace(L0,L+L0,N)    
y=np.linspace(L0,L+L0,N)

fig=plt.figure()
ax=fig.add_subplot(111)

color=2*np.log(np.hypot(Ex, Ey))
ax.streamplot(x,y,Ex,Ey,color=color,linewidth=1, cmap=plt.cm.inferno,density=2,arrowstyle='->',arrowsize=1.5)
plt.title('Campo eléctrico obtenido resolviendo la ecuación de Poisson')

for q, pos in cargas:
    ax.add_artist(Circle(pos, 0.05, color=charge_colors[q>0]))

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.set_aspect('equal')

X,Y=np.meshgrid(x,y)

fig = plt.figure()
plt.plot(y,V_temp)
plt.xlabel('$x$')
plt.ylabel('$V$') 
plt.title('Potencial a lo largo del eje de abcisas')

fig = plt.figure()
plt.plot(x,np.transpose(V_temp))
plt.xlabel('$y$')  #ax.set_xlabel('$x$')
plt.ylabel('$V$') #ax.set_ylabel('$rho$')
plt.title('Potencial a lo largo del eje de ordenadas')

fig = plt.figure()
erroresG=[0.0001,0.001,0.01,0.1,1]
iteracionesG=[1951,470,48,6,1]
plt.plot(np.log(erroresG),np.log(iteracionesG),'*')
plt.xlabel('Errores')
plt.ylabel('Iteraciones')

fig = plt.figure()
tamañoG=[5,7,10,20,30,40,50,60,70,80]
iteracionesG=[8,16,68,199,323,225,450,470,476,361]
plt.plot(tamañoG,iteracionesG,'*')
plt.xlabel('Tamaño (NxN)')
plt.ylabel('Iteraciones')