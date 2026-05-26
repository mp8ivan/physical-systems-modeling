import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

Nx = int(input("Numero de celdas unidad en la direccion x: "))
Ny = int(input("Numero de celdas unidad en la direccion y: "))
Nz = int(input("Numero de celdas unidad en la direccion z: "))

#Definimos los vectores primitivos:
def vectores_diamante(a):
    a1 = np.array([a,0,0])
    a2 = np.array([0,a,0])
    a3 = np.array([0,0,a])
    return np.column_stack((a1,a2,a3)) #esto nos devolvera una matriz formada por los tres vectores

A = vectores_diamante(1)

base_diamante = [np.array([0,0,0]),np.array([0.5,0.5,0]),np.array([0.5,0,0.5]),np.array([0,0.5,0.5]),
np.array([0.25,0.25,0.25]),np.array([0.75,0.75,0.25]),np.array([0.75,0.25,0.75]),np.array([0.25,0.75,0.75])]

#Creacion de la red:
#1. Formamos vectores n = [nx,ny,nz].
#2. Cada vector formado se multiplica por A.
#3. Añadimos el resultado a una lista.
red = []
for nx in range(0,Nx,1):
    for ny in range(0,Ny,1):
        for nz in range(0,Nz,1):
            n = np.array([nx,ny,nz])   
            for b in base_diamante:
                red.append(A.dot(n+b))
        
red = np.array(red)

print('Numero total de posiciones generadas: {}'.format(len(red)))

#Guardaremos los resultados en un fichero:
fichero = open('redCubica.xyz','w')
fichero.write('{}\n'.format(len(red)))
fichero.write('Red b.c.c.'+'\n')

for posicion in red:
    fichero.write('X'+' '+str(posicion[0])+' '+str(posicion[1])+' '+str(posicion[2])+' '+'\n') #X como elemento generico
fichero.close()

#Creamos la figura y los ejes 3D:
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Dibujamos los puntos:
ax.scatter(red[:,0],red[:,1],red[:,2],c='blue',s=50) #red[:,0] son todas nuestras coordenadas x, red[:,1] son todas nuestras coordenadas y y red[:,2] son todas nuestras coordenadas z
        
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_box_aspect([1, 1, 1]) #con esto, los ejes tendran la misma escala

plt.tight_layout()
plt.show()