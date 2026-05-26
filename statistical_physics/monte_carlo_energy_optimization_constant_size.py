import numpy as np
import matplotlib.pyplot as plt

N = int(input("Numero de particulas: "))
kT = [10,20,30,40]
M = int(input("Numero total de pasos del metodo Monte Carlo: "))

def E(nx,ny,nz):
    energia = nx**2 + ny**2 + nz**2 #simplificaremos calculos considerando la constante = 1
    return energia

for temperatura in kT: #realizamos la simulacion para cada uno de los valores de kT

    particulas = []
    for i in range(N):
        particulas.append([1,1,1])
    
    E_total = E(1,1,1)*N
    energias = [E_total]
    
    for paso in range(M):
    
        i = np.random.randint(0,N) #escogemos una particula aleatoria
        j = np.random.randint(0,3) #escogemos una componente aleatoria
        
        cambio = np.random.choice([-1,1])
        if particulas[i][j] + cambio >= 1:
            energia_vieja = E(*particulas[i]) #energia inicial de la particula a estudiar. E(*particulas[i]) equivale a E(particulas[0],particulas[1],...,particulas[i]))
            
            particula_propuesta = particulas[i].copy() #realizamos una copia de la i-esima particula
            particula_propuesta[j] = particulas[i][j] + cambio #aplicamos el cambio sobre la j-esima componente de dicha particula
            energia_nueva = E(*particula_propuesta) #nueva energia de la particula propuesta
        
            deltaE = energia_nueva - energia_vieja
            
            if deltaE <= 0 or np.random.rand() < np.exp(-deltaE/temperatura): #criterio de Metropolis
                particulas[i][j] = particulas[i][j] + cambio
                E_total += deltaE
    
        energias.append(E_total)
    plt.plot(energias, label='$kT = {}$'.format(temperatura))
    
#Evolucion de la energia total del sistema en funcion del numero de pasos de Monte Carlo:    
plt.xlabel("Pasos de Monte Carlo")
plt.ylabel("Energia total del sistema")
plt.legend()
plt.grid()
plt.show()