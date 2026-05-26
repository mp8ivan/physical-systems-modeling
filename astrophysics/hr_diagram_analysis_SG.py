import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit

#Funcion de normalizacion
def normalizar_espectro(wave, flux, iteraciones=5, grado_polinomio=4):
    mask_zoom = (wave >= 3900) & (wave <= 5000)
    w_zoom = wave[mask_zoom]
    f_zoom = flux[mask_zoom]
    w_fit = w_zoom.copy()
    f_fit = f_zoom.copy()
    
    for _ in range(iteraciones):
        coefs = np.polyfit(w_fit, f_fit, grado_polinomio)
        continuo_temp = np.polyval(coefs, w_zoom)
        
        residuos = f_zoom - continuo_temp
        desviacion = np.std(residuos)
  
        mask_buenos = (f_zoom - continuo_temp) > (-0.5 * desviacion)
        
        w_fit = w_zoom[mask_buenos]
        f_fit = f_zoom[mask_buenos]
        
    coefs_finales = np.polyfit(w_fit, f_fit, grado_polinomio)
    modelo_continuo_final = np.polyval(coefs_finales, w_zoom)
    
    return w_zoom, f_zoom / modelo_continuo_final


#Introducimos la estrella del catalogo a analizar:
estrella_id = input("Introduce el nombre del archivo MILES: ")
archivo_fits = f"miles_espectros/{estrella_id}.fits"

#Extraemos los datos del archivo FITS:
hdul = fits.open(archivo_fits)
data = hdul[0].data.flatten() 
header = hdul[0].header
hdul.close()

crval1 = header['CRVAL1']
cdelt1 = header['CDELT1']
naxis1 = header['NAXIS1']
wave = crval1 + cdelt1 * np.arange(naxis1)

flux = data

#Aplicamos la normalizacion al espectro:
w_norm, f_norm = normalizar_espectro(wave, flux)

#Espectro de la estrella (ampliado a las lineas de interes):
plt.figure(figsize=(12, 6))
plt.plot(w_norm, f_norm, color='black', linewidth=1, label=f'Espectro {estrella_id}')
plt.axhline(y=1.0, color='blue', linestyle=':', label='Continuo teórico (1.0)')
plt.axvline(x=4340.5, color='red', linestyle='--', alpha=0.8, label=r'H$\gamma$ (4340 Å)')
plt.axvline(x=4471.5, color='orange', linestyle='--', alpha=0.8, label='He I (4471 Å)')
plt.axvline(x=4481.2, color='cyan', linestyle='--', alpha=0.8, label='Mg II (4481 Å)')
plt.title(f"Espectro normalizado de la estrella")
plt.xlabel("Longitud de onda (Å)")
plt.ylabel("Flujo normalizado")
plt.xlim(4300, 4500) 
plt.ylim(0.4, 1.1)
plt.legend(loc="lower right")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

#Perfil gaussiano:
def gauss(x, A, mu, sigma):
    return 1.0 - A * np.exp(-0.5 * ((x - mu) / sigma)**2)

#-----------------LINEA DE BALMER-----------------

#Recortamos los datos alrededor de la linea de hidrogeno:
mascara_h = (w_norm >= 4330) & (w_norm <= 4350)
w_h = w_norm[mascara_h]
f_h = f_norm[mascara_h]

#Estimacion inicial a partir de la grafica:
p0_h = [0.35, 4340.5, 1.5]

#Ajuste con curve-fit:
popt_h, pcov_h = curve_fit(gauss, w_h, f_h, p0=p0_h)
A_h_opt, mu_h_opt, sigma_h_opt = popt_h

#Generamos la curva suave:
flujo_ajustado_h = gauss(w_h, *popt_h)

#Calculamos FWHM y EW:
fwhm_h = 2.355 * abs(sigma_h_opt)
ew_h = np.trapz(1.0 - flujo_ajustado_h, w_h)

resultados_h = f"Centro ajustado: {mu_h_opt:.2f} Å\nProfundidad: {A_h_opt:.3f}\nFWHM: {fwhm_h:.3f} Å\nEW: {ew_h:.3f} Å"

#Visualizacion del ajuste:
plt.figure(figsize=(8, 5))
plt.plot(w_h, f_h, 'o', color='black', markersize=3, label='Datos observados')
plt.plot(w_h, flujo_ajustado_h, '-', color='purple', linewidth=2.5, alpha=0.8, label='Ajuste gaussiano')
plt.axhline(y=1.0, color='blue', linestyle=':', label='Continuo (1.0)')
plt.title(r'Ajuste de la línea de H$\gamma$ (4340.5 $\AA$)')
plt.xlabel("Longitud de onda (Å)")
plt.ylabel("Flujo normalizado")
plt.legend()
plt.text(0.95, 0.05, resultados_h, transform=plt.gca().transAxes, 
         fontsize=10, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='lightgray', alpha=0.8))
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

#-----------------LINEA DE HE I-----------------

#Recortamos los datos alrededor del helio neutro:
mascara_he = (w_norm >= 4465) & (w_norm <= 4478)
w_he = w_norm[mascara_he]
f_he = f_norm[mascara_he]

#Estimacion inicial a partir de la grafica:
p0_he = [0.17, 4471.5, 1.0]

#Ajuste con curve_fit:
popt_he, pcov_he = curve_fit(gauss, w_he, f_he, p0=p0_he)
A_opt, mu_opt, sigma_opt = popt_he

#Generamos la curva suave:
flujo_ajustado_he = gauss(w_he, *popt_he)

#Calculamos FWHM y EW:
fwhm_he = 2.355 * abs(sigma_opt) 
ew_he = np.trapz(1.0 - flujo_ajustado_he, w_he)

resultados_he = f"Centro ajustado: {mu_opt:.2f} Å\nProfundidad: {A_opt:.3f}\nFWHM: {fwhm_he:.3f} Å\nEW: {ew_he:.3f} Å"

#Visualizacion del ajuste:
plt.figure(figsize=(8, 5))
plt.plot(w_he, f_he, 'o', color='black', markersize=4, label='Datos del catálogo')
plt.plot(w_he, flujo_ajustado_he, '-', color='red', linewidth=2.5, alpha=0.8, label='Ajuste gaussiano')
plt.axhline(y=1.0, color='blue', linestyle=':', label='Continuo (1.0)')
plt.title(r'Ajuste de la línea de He I (4471.5 $\AA$)')
plt.xlabel("Longitud de onda (Å)")
plt.ylabel("Flujo normalizado")
plt.legend()
plt.text(0.95, 0.05, resultados_he, transform=plt.gca().transAxes, 
         fontsize=10, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='lightgray', alpha=0.8))
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

#-----------------LINEA DE MG II-----------------

#Recortamos los datos alrededor de la linea de magnesio ionizado:
mascara_mg = (w_norm >= 4478) & (w_norm <= 4485)
w_mg = w_norm[mascara_mg]
f_mg = f_norm[mascara_mg]

#Estimacion inicial a partir de la grafica:
p0_mg = [0.13, 4481.2, 0.8]

#Ajuste con curve-fit:
popt_mg, pcov_mg = curve_fit(gauss, w_mg, f_mg, p0=p0_mg)
A_mg_opt, mu_mg_opt, sigma_mg_opt = popt_mg

#Generamos la curva suave:
flujo_ajustado_mg = gauss(w_mg, *popt_mg)

#Calculamos FWHM y EW:
fwhm_mg = 2.355 * abs(sigma_mg_opt)
ew_mg = np.trapz(1.0 - flujo_ajustado_mg, w_mg)

resultados_mg = f"Centro ajustado: {mu_mg_opt:.2f} Å\nProfundidad: {A_mg_opt:.3f}\nFWHM: {fwhm_mg:.3f} Å\nEW: {ew_mg:.3f} Å"

#Visualizacion del ajuste:
plt.figure(figsize=(8, 5))
plt.plot(w_mg, f_mg, 'o', color='black', markersize=4, label='Datos observados')
plt.plot(w_mg, flujo_ajustado_mg, '-', color='cyan', linewidth=2.5, alpha=0.8, label='Ajuste gaussiano')
plt.axhline(y=1.0, color='blue', linestyle=':', label='Continuo (1.0)')
plt.title(r'Ajuste de la línea de Mg II (4481.2 $\AA$)')
plt.xlabel("Longitud de onda (Å)")
plt.ylabel("Flujo normalizado")
plt.legend()
plt.text(0.95, 0.05, resultados_mg, transform=plt.gca().transAxes, 
         fontsize=10, verticalalignment='bottom', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='lightgray', alpha=0.8))
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()