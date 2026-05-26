import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#Definimos las constantes y los parametros fisicos:
X = 0.7
Y = 0.28
Z = 0.02
gamma = 5/3
k_bol = 1.380649e-16 
mh = 1.6726e-24      
mu = 1/(2*X+0.75*Y+0.5*Z)

#Constantes del Sol (C.G.S.):
M_sol = 1.989e33
R_sol = 6.957e10
G_cte = 6.674e-8
L_sol = 3.828e33

#Parametros combinados para simplificar ecuaciones:
C = (G_cte*M_sol*mh)/(1e6*R_sol*k_bol)
D = 1/C

#Definimos la funcion de generacion de energia:
def epsilon(rho,T6):
    """
    rho: float. Densidad actual.
    T6: float. Temperatura actual (en millones de K).
    La funcion devolvera la tasa de generacion de energia combinada (pp + CNO).
    """
    ePP_tilde = 9.3e6 * (X**2) * rho * (T6**(-2/3)) * np.exp(-33.8 / (T6**(1/3)))
    eCNO_tilde = 2.5e28 * X * Z * rho * (T6**(-2/3)) * np.exp(-152.3 / (T6**(1/3)))
    return ePP_tilde + eCNO_tilde

#Definimos el sistema de ecuaciones diferenciales:
def Eddington_Modificado(x,Y):
    """
    x: float. Variable independiente (ln(P_tilde)).
    Y: array. Vector con las variables dependientes [r_tilde,M_r_tilde,L_r_tilde,ln(T6)].
    """
    r = Y[0]
    M_r = Y[1]
    L_r = Y[2]
    lnT6 = Y[3]
    
    P_tilde = np.exp(x)
    T6 = np.exp(lnT6)
    rho = (P_tilde*mu)/(D*T6)
    
    k_tilde = 1e6*rho*(T6**(-3.5))
    e_tilde = epsilon(rho,T6)
    
    n_rad = (k_tilde*L_r*P_tilde)/(4*np.pi*M_r*(T6**4))
    n_conv = 1-(1/gamma)
    n = min(n_rad, n_conv)
    
    dr_dx = -((r**2)*P_tilde)/(M_r*rho)
    dM_dx = -(4*np.pi*(r**4)*P_tilde)/M_r
    dL_dx = -(4*np.pi*(r**4)*P_tilde*e_tilde)/M_r
    dlnT6_dx = n
    
    return [dr_dx,dM_dx,dL_dx,dlnT6_dx]

#Proceso de resolucion:
P0 = 17.28
T60 = 16.1

r_c = 1e-4
rho_c = (P0*mu)/(D*T60)
M_c = (4/3)*np.pi*rho_c*(r_c**3)
L_c = M_c * epsilon(rho_c, T60)

Y0 = [r_c,M_c,L_c,np.log(T60)]
x_inicial = np.log(P0)
Ps = P0*1e-8 
x_final = np.log(Ps)

solucion = solve_ivp(Eddington_Modificado,t_span=[x_inicial, x_final],y0=Y0,method='RK45',rtol=1e-8,atol=1e-8)

ln_P = solucion.t
P_res = np.exp(ln_P)
r_res = solucion.y[0]
M_res = solucion.y[1]
L_res = solucion.y[2]
T6_res = np.exp(solucion.y[3])
rho_res = (P_res*mu)/(D*T6_res)

#Conversion a C.G.S.:
Unidad_P = (G_cte*(M_sol**2))/(R_sol**4)
r_CGS = r_res*R_sol
M_CGS = M_res*M_sol
L_CGS = L_res*L_sol
P_CGS = P_res*Unidad_P
T_CGS = T6_res*1e6

#Cálculo de gradientes y ratio:
k_vec = 1e6*rho_res*(T6_res**(-3.5))
nabla_rad = (k_vec*L_res*P_res)/(4*np.pi*M_res*(T6_res**4))
nabla_conv = np.full_like(r_res, 1-(1/gamma))
factor_conversion = T6_res*(M_res*rho_res)/((r_res**2)*P_res)
dTdr_rad = factor_conversion*nabla_rad
dTdr_conv = factor_conversion*nabla_conv

#Primer bloque de graficas:
for data,label in zip([r_res,M_res,L_res,T6_res],['$\\tilde{r}(R_{\\odot})$','$\\tilde{M}_r(M_{\\odot})$','$\\tilde{L}_r(L_{\\odot})$','$T_6$']):
    plt.figure()
    plt.plot(ln_P,data,color='blue',linewidth=2.5)
    plt.xlabel('$\\ln(\\tilde{P})$',fontsize=11) 
    plt.ylabel(label,fontsize=11)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

#Segundo bloque de graficas:
for data,label in zip([M_res, P_res, L_res, T6_res],['$\\tilde{M}_r(M_{\\odot})$','$\\tilde{P}$','$\\tilde{L}_r(L_{\\odot})$','$T_6$']):
    plt.figure()
    plt.plot(r_res,data,color='blue',linewidth=2.5)
    plt.xlabel('$\\tilde{r}(R_{\\odot})$',fontsize=11)
    plt.ylabel(label,fontsize=11)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

#Tercer bloque de graficas:
for data,label,unit in zip([M_CGS,P_CGS,L_CGS,T_CGS],['M_r','P','L_r','T'],['g','dinas/cm$^2$','erg/s','K']):
    plt.figure()
    plt.plot(r_CGS,data,color='blue',linewidth=2.5)
    plt.xlabel('$r$ (cm)',fontsize=11) 
    plt.ylabel(f'${label}$ ({unit})',fontsize=11)
    plt.grid(True,linestyle='--',alpha=0.7)
    plt.tight_layout()

#Criterio de Schwarzschild:
delta_nabla = nabla_rad - nabla_conv

plt.figure()
plt.plot(r_res, nabla_rad, color='blue', linewidth=2.5, label='$\\nabla_{rad}$ (Radiativo)')
plt.plot(r_res, nabla_conv, color='red', linestyle='--', linewidth=2.0, label='$\\nabla_{ad}$ (Adiabatico)')
plt.xlabel('$\\tilde{r}(R_{\\odot})$', fontsize=11)
plt.ylabel('$\\nabla$', fontsize=11)
plt.xlim(0, 0.4) 
plt.ylim(0, 1.5)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(r_res, delta_nabla, color='purple', linewidth=2.5, label='$\\Delta \\nabla = \\nabla_{rad} - \\nabla_{ad}$')
plt.axhline(y=0.0, color='black', linestyle='-', linewidth=2.0, label='Limite de estabilidad (0)')
plt.fill_between(r_res, 0, delta_nabla, where=(delta_nabla <= 0), color='skyblue', alpha=0.3, label='Zona estable (radiacion)')
plt.fill_between(r_res, 0, delta_nabla, where=(delta_nabla > 0), color='orange', alpha=0.3, label='Zona inestable (conveccion)')
plt.xlabel('$\\tilde{r}(R_{\\odot})$', fontsize=11)
plt.ylabel('$\\Delta \\nabla$', fontsize=11)
plt.xlim(0, 0.4)
plt.ylim(-0.5, 2.0)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()

plt.show()