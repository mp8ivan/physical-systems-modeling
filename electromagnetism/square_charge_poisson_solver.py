import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
charge_colors = {True: '#aa0000', False: '#0000aa'}

N=61
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

rho = np.zeros((N, N))

lado_cuadrado = 1.5
M = 500
carga_por_punto = 3.

for k in range(M):
    frac = k / M
    if 0 <= frac < 0.25:
        x = -lado_cuadrado + 2 * lado_cuadrado * frac / 0.25
        y = -lado_cuadrado
    elif 0.25 <= frac < 0.5:
        x = lado_cuadrado
        y = -lado_cuadrado + 2 * lado_cuadrado * (frac - 0.25) / 0.25
    elif 0.5 <= frac < 0.75:
        x = lado_cuadrado - 2 * lado_cuadrado * (frac - 0.5) / 0.25
        y = lado_cuadrado
    else: 
        x = -lado_cuadrado
        y = lado_cuadrado - 2 * lado_cuadrado * (frac - 0.75) / 0.25
    i, j = Pos((x, y))
    if 0 <= i < N and 0 <= j < N:
        rho[i, j] += carga_por_punto

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

ax.add_patch(Rectangle((-lado_cuadrado, -lado_cuadrado), 2 * lado_cuadrado, 2 * lado_cuadrado, edgecolor=charge_colors[True], 
                       facecolor='none', linewidth=2))

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
ax.set_aspect('equal')

X,Y=np.meshgrid(x,y)

plt.show()