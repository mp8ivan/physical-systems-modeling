import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Definimos las constantes y los parametros fisicos:
G = 6.67430e-8
c = 2.99792458e10
h = 6.62607015e-27
hbar = h/(2*np.pi)
m_e = 9.10938356e-28
m_B = 1.6726e-24
M_sol = 1.98847e33
R_tierra = 6.371e8
mu_e = 2
lambdabar = hbar/(m_e*c)

#Momento de Fermi:
def x_F(rho):
    """
    rho: float. Densidad actual.
    La funcion devolvera el momento de Fermi correspondiente a esa densidad.
    """
    return (rho*3*(lambdabar**3)*(np.pi**2)/(mu_e*m_B))**(1/3)

#Funcion de Chandrasekhar:
def phi(x):
    """
    x: float. Momento adimensional.
    La funcion devolvera el valor de la funcion de Chandrasekhar para ese x.
    """
    numerador = x*np.sqrt(1+x**2)*((2*x**2)/3-1)+np.log(x+np.sqrt(1+x**2))
    return numerador/(8*np.pi**2)

#Presion degenerada:
def P(x):
    """
    x: float. Momento adimensional.
    La funcion devolvera la presion degenerada correspondiente a ese x.
    """
    return (m_e*c**2)*phi(x)/(lambdabar**3)

#Indice adiabatico local:
def gamma(x):
    """
    x: float. Momento adimensional.
    La funcion devolvera el indice adiabatico local correspondiente a ese x.
    """
    denominador = 9*np.pi**2*np.sqrt(1+x**2)*phi(x)
    return (x**5)/denominador

#Definimos el sistema de ecuaciones diferenciales:
def enanaBlanca(x,Y):
    """
    x: float. Variable independiente (ln(rho_tilde)).
    Y: array. Vector con las variables dependientes [r_tilde,M_r_tilde].
    La función devolverá dY_dx: list. Derivadas de las variables dependientes con respecto a x.
    """
    r_tilde = Y[0]
    Mr_tilde = Y[1]
    
    #1. Obtenemos la densidad normalizada:
    rho_tilde = np.exp(x)

    #2. Calculamos la presion en C.G.S. para obtener x_F:
    rho = rho_tilde*(M_sol/(R_tierra**3))

    #3. Calculamos la presion en C.G.S. y el indice adiabatico local:
    P_cgs = P(x_F(rho))
    gamma_cgs = gamma(x_F(rho))

    #4. Normalizamos la presion:
    P_tilde = P_cgs/(G*(M_sol**2)/(R_tierra**4))

    #5. Calculamos las derivadas:
    dr_tilde_dx = -((r_tilde**2)*P_tilde*gamma_cgs)/(Mr_tilde*rho_tilde)
    dM_tilde_dx = -(4*np.pi*(r_tilde**4)*P_tilde*gamma_cgs)/(Mr_tilde)
    
    return [dr_tilde_dx,dM_tilde_dx]

#Condiciones iniciales:
rho_centro = 1e6 #densidad central de prueba
rho_superficie = 1
r_centro = 1e-6*R_tierra
Mr_centro = (4/3)*np.pi*(r_centro**3)*rho_centro

#Normalizamos las condiciones iniciales:
rho_tilde_centro = rho_centro/(M_sol/(R_tierra**3))
rho_tilde_superficie = rho_superficie/(M_sol/(R_tierra**3))
r_tilde_centro = r_centro/R_tierra
M_tilde_centro = Mr_centro/M_sol

#Resolvemos el sistema a partir de Runge-Kutta:
sol = solve_ivp(enanaBlanca,t_span=[np.log(rho_tilde_centro),np.log(rho_tilde_superficie)],y0=[r_tilde_centro,M_tilde_centro],method='RK45',rtol=1e-8,atol=1e-8)

#Extraemos los resultados:
ln_rho_tilde = sol.t
rho_tilde = np.exp(ln_rho_tilde)
r_tilde = sol.y[0]
Mr_tilde = sol.y[1]

#Conversion a C.G.S.:
rho_cgs = rho_tilde*(M_sol/(R_tierra**3))
r_cgs = r_tilde*R_tierra
Mr_cgs = Mr_tilde*M_sol
P_cgs = P(x_F(rho_cgs))
gamma_vec = gamma(x_F(rho_cgs))

#Imprimimos los resultados:
print(f"Densidad central inicial: {rho_centro} g/cm^3")
print(f"Masa final: {Mr_cgs[-1]:.2e} g")
print(f"Radio final: {r_cgs[-1]:.2e} cm")

#Primer bloque de graficas:
for data,label in zip([r_tilde,Mr_tilde],['$\\tilde{r} (R_{\\oplus})$','$\\tilde{M}_r (M_{\\odot})$']):
    plt.figure()
    plt.plot(ln_rho_tilde,data,color='blue',linewidth=2.5)
    plt.xlabel('$\\ln\\left(\\frac{\\rho}{M_{\\odot}/R_{\\oplus}^3}\\right)$',fontsize=11) 
    plt.ylabel(label,fontsize=11)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

#Segundo bloque de graficas:
for data,label in zip([Mr_tilde,P_cgs/(G*(M_sol**2)/(R_tierra**4))],['$\\tilde{M}_r (M_{\\odot})$','$\\frac{P}{M_{\\odot}/R_{\\oplus}^3}$']):
    plt.figure()
    plt.plot(r_tilde,data,color='blue',linewidth=2.5)
    plt.xlabel('$\\tilde{r} (R_{\\oplus})$',fontsize=11)
    plt.ylabel(label,fontsize=11)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

#Tercer bloque de graficas:
for data,label,unit in zip([Mr_cgs,P_cgs],['M_r','P'],['g','dyn/cm$^2$']):
    plt.figure()
    plt.plot(r_cgs,data,color='blue',linewidth=2.5)
    plt.xlabel('$r$ (cm)',fontsize=11) 
    plt.ylabel(f'${label}$ ({unit})',fontsize=11)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

#Cuarto bloque de graficas:
plt.figure()
plt.plot(r_cgs,rho_cgs,color='blue',linewidth=2.5)
plt.xlabel('$r$ (cm)',fontsize=11)
plt.ylabel('$\\rho$ (g/cm$^3$)',fontsize=11)
plt.grid(True,linestyle='--',alpha=0.7)
plt.tight_layout()

plt.figure()
plt.plot(x_F(rho_cgs),gamma_vec,color='blue',linewidth=2.5)
plt.xlabel('$x_F$',fontsize=11)
plt.ylabel('$\\gamma$',fontsize=11)
plt.grid(True,linestyle='--',alpha=0.7)
plt.tight_layout()

plt.figure()
plt.plot(rho_cgs,P_cgs,color='blue',linewidth=2.5)
plt.xlabel('$\\rho$ (g/cm$^3$)',fontsize=11)
plt.ylabel('$P$ (dyn/cm$^2$)',fontsize=11)
plt.grid(True,linestyle='--',alpha=0.7)
plt.tight_layout()

#Quinto bloque de graficas:
plt.figure()
plt.plot(r_cgs,Mr_cgs,color='blue',linewidth=2.5)
plt.xlabel('$r$ (cm)',fontsize=11)
plt.ylabel('$M_r$ (g)',fontsize=11)
plt.grid(True,linestyle='--',alpha=0.7)
plt.tight_layout()

densidades_centrales = np.logspace(4,15,40)
radios = []
masas = []
for densidad_centro in densidades_centrales:
    #Condiciones iniciales:
    densidad_superficie = 1
    Masa_centro = (4/3)*np.pi*(r_centro**3)*densidad_centro

    #Normalizamos las condiciones iniciales:
    densidad_tilde_centro = densidad_centro/(M_sol/(R_tierra**3))
    densidad_tilde_superficie = densidad_superficie/(M_sol/(R_tierra**3))
    Masa_tilde_centro = Masa_centro/M_sol

    #Resolvemos el sistema a partir de Runge-Kutta:
    solucion = solve_ivp(enanaBlanca,t_span=[np.log(densidad_tilde_centro),np.log(densidad_tilde_superficie)],y0=[r_tilde_centro,Masa_tilde_centro],method='RK45',rtol=1e-8,atol=1e-8)

    #Extraemos la solucion:
    radios.append(solucion.y[0][-1])
    masas.append(solucion.y[1][-1])

#Conversion a C.G.S.:
radios_cgs = np.array(radios)*R_tierra

plt.figure()
plt.plot(radios_cgs,masas,color='blue',linewidth=2.5)
plt.xlabel('$r$ (cm)',fontsize=11)
plt.ylabel('$M_r$ (g)',fontsize=11)
plt.grid(True,linestyle='--',alpha=0.7)
plt.tight_layout()

plt.figure()
plt.plot(densidades_centrales,masas,color='blue',linewidth=2.5)
plt.xlabel('$\\rho_0$ (g/cm$^3$)',fontsize=11)
plt.ylabel('$M_r$ (g)',fontsize=11)
plt.axhline(1.44,linestyle='--',color='red',label='Limite de Chandrasekhar')
plt.xscale('log')
plt.grid(True,linestyle='--',alpha=0.7)
plt.legend()
plt.tight_layout()