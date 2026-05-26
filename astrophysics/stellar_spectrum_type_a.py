import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import splrep, splev

#Funcion de normalizacion:
def normalizar_espectro(wave, flux):
    mask_zoom = (wave >= 3900) & (wave <= 5000)
    w_zoom = wave[mask_zoom]
    f_zoom = flux[mask_zoom]
    
    zonas_continuo = [
    (4040, 4070),  
    (4150, 4190),  
    (4250, 4290),  
    (4420, 4460),  
    (4530, 4580),  
    (4740, 4790),  
    (4950, 5000)   
]
    
    w_puntos, f_puntos = [], []
    for zona in zonas_continuo:
        mask_zona = (w_zoom >= zona[0]) & (w_zoom <= zona[1])
        if np.any(mask_zona):
            w_puntos.append(np.median(w_zoom[mask_zona]))
            f_puntos.append(np.percentile(f_zoom[mask_zona], 95))
            
    spline = splrep(w_puntos, f_puntos, k=1, s=0)
    modelo_continuo = splev(w_zoom, spline)
    return w_zoom, f_zoom / modelo_continuo

#Carga de datos:
wave_prob, flux_prob = np.loadtxt("estrella3.dat", skiprows=3, unpack=True)
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
plt.plot(w_prob_norm, f_prob_norm, color='black', label='Estrella 3')
plt.axvline(x=3933.7, color='green', linestyle='--', label='Ca II K (3933 Å)')
plt.axvline(x=3968.5, color='green', linestyle='--', label='Ca II H (3969 Å)')
plt.axvline(x=4128.0, color='gray', linestyle='--', label='Si II (4128 Å)')
plt.axvline(x=4130.9, color='gray', linestyle='--', label='Si II (4131 Å)')
plt.axvline(x=4233.2, color='brown', linestyle='--', label='Fe II (4233 Å)')
plt.axvline(x=4481.2, color='cyan', linestyle='--', label='Mg II (4481 Å)')
plt.axvline(x=4703.0, color='magenta', linestyle='--', label='Mg I (4703 Å)')
plt.title(f"Comparación espectral: estrella 3 vs {tipo_f}")
plt.legend(loc="lower right")
plt.xlim(3900, 5010); plt.ylim(0.0, 1.15)
plt.tight_layout()
plt.show()

#Segunda grafia (lineas de Balmer):
plt.figure(figsize=(12, 6))
plt.plot(w_m_norm, f_m_norm, color='red', alpha=0.6, label=f'MILES ({tipo_f})')
plt.plot(w_prob_norm, f_prob_norm, color='black', label='Estrella 3')

balmer = {
    r'$H\epsilon$ (3970.1 Å)': 3970.1,
    r'$H\delta$ (4101.7 Å)': 4101.7,
    r'$H\gamma$ (4340.5 Å)': 4340.5,
    r'$H\beta$ (4861.3 Å)': 4861.3
}

colors = ['purple', 'blue', 'green', 'darkred']
for i, (name, wave) in enumerate(balmer.items()):
    plt.axvline(x=wave, color=colors[i], linestyle='--', alpha=0.8, label=name)

plt.title(f"Comparación espectral: estrella 3 vs {tipo_f}")
plt.xlabel("Longitud de onda (Å)"); plt.ylabel("Flujo normalizado")
plt.legend(loc="lower right", ncol=3, fontsize='small') 
plt.xlim(3900, 5010); plt.ylim(0.0, 1.15)
plt.tight_layout()
plt.show()