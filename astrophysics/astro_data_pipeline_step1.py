import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Definimos las constantes y los parametros fisicos:
mu = 0.62
beta = 0.9995
G = 6.67430e-8
M_sol = 1.98847e33
R_sol = 6.957e10

#Calculamos la constante politropica adimensional:
K_tilde = 2.5*((1-beta)/((mu**4)*(beta**4)))**(1/3)

#Definimos el sistema de ecuaciones diferenciales:
def Eddington(x,Y):

    """
    x: float. Variable independiente (ln(P_tilde)).
    Y: array. Vector con las variables dependientes [r_tilde,M_tilde].
    La función devolverá dY_dx: list. Derivadas de las variables dependientes con respecto a x.
    """

    #1. Desempaquetamos las variables dependientes:
    r = Y[0]
    M_r = Y[1]
    
    #2. Calculamos la densidad y la presion actuales a partir del cambio de variable:
    P_tilde = np.exp(x)
    rho_tilde = (P_tilde/K_tilde)**(3/4)
    
    #3. Calculamos las derivadas de r_tilde y M_r_tilde con respecto a x:
    dr_dx = -((r**2)*P_tilde)/(M_r*rho_tilde) #gradiente del radio
    dM_dx = -(4*np.pi*(r**4)*P_tilde)/M_r #gradiente de la masa
    
    #4. Devolvemos las derivadas como una lista:
    return [dr_dx,dM_dx]

#Bucle de busqueda a partir del metodo de la biseccion:
Pc_min = 1.0 #limite inferior de prueba
Pc_max = 1e5 #limite superior de prueba
tol = 1e-4 #queremos que el radio final difiera de 1 en menos de 0.0001
max_iter = 50 #maximo numero de iteraciones para evitar bucles infinitos

exito = False
for i in range(max_iter):
    
    #1. Proponemos como presion central el punto medio entre Pc_min y Pc_max:
    Pc_prueba = (Pc_min+Pc_max)/2
    
    #2. Calculamos las condiciones iniciales:
    r_c = 10**(-4)
    rho_c = (Pc_prueba/K_tilde)**(3/4)
    M_c = (4/3)*np.pi*(r_c**3)*rho_c
    
    #3. Definimos el vector de estado inicial:
    Y0 = [r_c,M_c]
    
    #4. Definimos el intervalo de integracion:
    x_inicial = np.log(Pc_prueba)
    Ps = 10**(-8)*Pc_prueba 
    x_final = np.log(Ps)
    
    #5. Resolvemos el sistema de ecuaciones diferenciales por Runge-Kutta:
    solucion = solve_ivp(Eddington,t_span=[x_inicial,x_final],y0=Y0,method='RK45',rtol=1e-8,atol=1e-8)
    
    #6. Extraemos los datos en la superficie:
    r_superficie = solucion.y[0][-1]
    M_superficie = solucion.y[1][-1]
    
    #7. Calculamos cuanto nos hemos desviado del radio ideal:
    error = r_superficie-1 
    
    #8. Comprobamos si hemos acertado:
    if abs(error) < tol:
        Pc_optimo = Pc_prueba
        exito = True
        break #convergencia alcanzada
    
    #9. Corregimos la desviacion:
    if error > 0:
        Pc_min = Pc_prueba #el radio es muy grande, necesitamos aumentar la presion central para compactar la estrella
    else:
        Pc_max = Pc_prueba #el radio es muy pequeño, necesitamos disminuir la presion central para expandir la estrella
else: 
    print("No se alcanzo la convergencia despues de {} iteraciones".format(max_iter))

#Visualizamos los resultados:
if exito == True:
    
    #1. Extraemos los datos de la solucion final:
    ln_P = solucion.t
    r_tilde = solucion.y[0]
    M_tilde = solucion.y[1]
    rho_tilde = (np.exp(ln_P)/K_tilde)**(3/4)

    #2. Primera grafica:
    plt.figure()
    plt.plot(ln_P,r_tilde,color='blue',linewidth=2.5)
    plt.xlabel('$\\ln(\\tilde{P})(M_{\odot},R_{\odot})$',fontsize=11)
    plt.ylabel('$\\tilde{r}(R_{\odot})$',fontsize=11)
    plt.title('Modelo de Eddington (unidades normalizadas)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #3. Segunda grafica:
    plt.figure() 
    plt.plot(ln_P,M_tilde,color='blue',linewidth=2.5)
    plt.xlabel('$\\ln(\\tilde{P})(M_{\odot},R_{\odot})$',fontsize=11)
    plt.ylabel('$\\tilde{M}_r(M_{\odot})$',fontsize=11)
    plt.title('Modelo de Eddington (unidades normalizadas)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #4. Tercera grafica:
    plt.figure() 
    plt.plot(r_tilde,M_tilde,color='blue',linewidth=2.5)
    plt.xlabel('$\\tilde{r}(R_{\odot})$',fontsize=11)
    plt.ylabel('$\\tilde{M}_r(M_{\odot})$',fontsize=11)
    plt.title('Modelo de Eddington (unidades normalizadas)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #5. Cuarta grafica:
    plt.figure() 
    plt.plot(r_tilde,np.exp(ln_P),color='blue',linewidth=2.5)
    plt.xlabel('$\\tilde{r}(R_{\odot})$',fontsize=11)
    plt.ylabel('$\\tilde{P}(M_{\odot},R_{\odot})$',fontsize=11)
    plt.title('Modelo de Eddington (unidades normalizadas)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #6. Quinta grafica:
    plt.figure() 
    plt.plot(r_tilde,rho_tilde,color='blue',linewidth=2.5)
    plt.xlabel('$\\tilde{r}(R_{\odot})$',fontsize=11)
    plt.ylabel('$\\tilde{\\rho}(M_{\odot},R_{\odot})$',fontsize=11)
    plt.title('Modelo de Eddington (unidades normalizadas)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #7. Sexta grafica:
    plt.figure() 
    plt.plot(rho_tilde,np.exp(ln_P),color='blue',linewidth=2.5)
    plt.xlabel('$\\tilde{\\rho}(M_{\odot},R_{\odot})$',fontsize=11)
    plt.ylabel('$\\tilde{P}(M_{\odot},R_{\odot})$',fontsize=11)
    plt.title('Modelo de Eddington (unidades normalizadas)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #8. Septima grafica:
    plt.figure() 
    plt.plot((M_sol/(R_sol**3))*rho_tilde,(G*(M_sol**2)/(R_sol**4))*np.exp(ln_P),color='blue',linewidth=2.5)
    plt.xlabel('$\\rho$ (g/cm$^3$)',fontsize=11)
    plt.ylabel('$P$ (dyn/cm$^2$)',fontsize=11)
    plt.title('Modelo de Eddington (C.G.S.)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #9. Octava grafica:
    plt.figure() 
    plt.plot(r_tilde*R_sol,M_tilde*M_sol,color='blue',linewidth=2.5)
    plt.xlabel('$r$ (cm)',fontsize=11)
    plt.ylabel('$M$ (g)',fontsize=11)
    plt.title('Modelo de Eddington (C.G.S.)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #10. Novena grafica:
    plt.figure() 
    plt.plot(r_tilde*R_sol,(G*(M_sol**2)/(R_sol**4))*np.exp(ln_P),color='blue',linewidth=2.5)
    plt.xlabel('$r$ (cm)',fontsize=11)
    plt.ylabel('$P$ (dyn/cm$^2$)',fontsize=11)
    plt.title('Modelo de Eddington (C.G.S.)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    #11. Decima grafica:
    plt.figure() 
    plt.plot(r_tilde*R_sol,(M_sol/(R_sol**3))*rho_tilde,color='blue',linewidth=2.5)
    plt.xlabel('$r$ (cm)',fontsize=11)
    plt.ylabel('$\\rho$ (g/cm$^3$)',fontsize=11)
    plt.title('Modelo de Eddington (C.G.S.)',fontsize=12)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

    plt.show()

    #Valores finales obtenidos para la estrella:
    print("Para P_c = {}:".format(Pc_optimo))
    print("---Unidades normalizadas---")
    print(f"Radio final (r~): {r_tilde[-1]:.6f}")
    print(f"Masa final (M~): {M_tilde[-1]:.6f}")

    print("---Unidades C.G.S.---")
    print(f"Radio final (r): {r_tilde[-1] * R_sol:.4e} cm")
    print(f"Masa final (M): {M_tilde[-1] * M_sol:.4e} g")

else:
    print("Error: no se encontro una solucion valida para graficar")