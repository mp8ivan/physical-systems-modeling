import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

#Funcion de normalizacion:
def normalizar_espectro(wave, flux, iteraciones=5, grado_polinomio=4):
    mask_zoom = (wave >= 3900) & (wave <= 5000)
    w_zoom = wave[mask_zoom]
    f_zoom = flux[mask_zoom]
    w_fit = w_zoom.copy()
    f_fit = f_zoom.copy()
    #Bucle de rechazo iterativo (Sigma-Clipping):
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

#Carga de datos:
wave_prob, flux_prob = np.loadtxt("estrella1.dat", skiprows=3, unpack=True)
w_prob_norm, f_prob_norm = normalizar_espectro(wave_prob, flux_prob)

#Inputs para que el usuario escoja la estrella del catalogo a comparar: 
nombre_f = input("Escribe el nombre del FITS a comparar: ")
tipo_f = input(f"¿Qué tipo espectral es '{nombre_f}'?: ")

#Cargamos y normalizamos el MILES:
with fits.open(f"miles_espectros/{nombre_f}.fits") as hdul:
    f_m = hdul[0].data.flatten()
    h = hdul[0].header
    w_m = h['CRVAL1'] + h['CDELT1'] * np.arange(h['NAXIS1'])

w_m_norm, f_m_norm = normalizar_espectro(w_m, f_m / np.median(f_m))

#Primera grafica (tipo espectral):
plt.figure(figsize=(12, 6))
plt.plot(w_m_norm, f_m_norm, color='red', alpha=0.6, label=f'MILES ({tipo_f})')
plt.plot(w_prob_norm, f_prob_norm, color='black', label='Estrella 1')
plt.axvline(x=3933.7, color='green', linestyle='--', label='Ca II K (3933 Å)')
plt.axvline(x=4026.2, color='orange', linestyle='--', label='He I (4026 Å)')
plt.axvline(x=4128.0, color='gray', linestyle='--', label='Si II (4128 Å)')
plt.axvline(x=4130.9, color='gray', linestyle='--', label='Si II (4131 Å)')
plt.axvline(x=4471.5, color='orange', linestyle='--', label='He I (4471 Å)')
plt.axvline(x=4481.2, color='cyan', linestyle='--', label='Mg II (4481 Å)')
plt.title(f"Comparación espectral: estrella 1 vs {tipo_f}")
plt.legend(loc="lower right")
plt.xlim(3900, 5010); plt.ylim(0.35, 1.15)
plt.tight_layout()
plt.show()

#Segunda grafia (lineas de Balmer):
plt.figure(figsize=(12, 6))
plt.plot(w_m_norm, f_m_norm, color='red', alpha=0.6, label=f'MILES ({tipo_f})')
plt.plot(w_prob_norm, f_prob_norm, color='black', label='Estrella 1')

balmer = {
    r'$H\epsilon$ (3970.1 Å)': 3970.1,
    r'$H\delta$ (4101.7 Å)': 4101.7,
    r'$H\gamma$ (4340.5 Å)': 4340.5,
    r'$H\beta$ (4861.3 Å)': 4861.3
}

colors = ['purple', 'blue', 'green', 'darkred']
for i, (name, wave) in enumerate(balmer.items()):
    plt.axvline(x=wave, color=colors[i], linestyle='--', alpha=0.8, label=name)

plt.title(f"Comparación espectral: estrella 1 vs {tipo_f}")
plt.xlabel("Longitud de onda (Å)"); plt.ylabel("Flujo normalizado")
plt.legend(loc="lower right", ncol=3, fontsize='small') 
plt.xlim(3900, 5010); plt.ylim(0.35, 1.15)
plt.tight_layout()
plt.show()