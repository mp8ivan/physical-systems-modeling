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
    (4000, 4020), 
    (4080, 4110), 
    (4220, 4240),  
    (4310, 4330),  
    (4440, 4470),  
    (4500, 4540),  
    (4600, 4650), 
    (4740, 4780),  
    (4880, 4920),  
    (4980, 5010)   
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
wave_prob, flux_prob = np.loadtxt("estrella4.dat", skiprows=3, unpack=True)
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
plt.plot(w_prob_norm, f_prob_norm, color='black', label='Estrella 4')
plt.axvline(x=3933.7, color='green', linestyle='--', label='Ca II K (3933 Å)')
plt.axvline(x=3968.5, color='green', linestyle='--', label='Ca II H (3969 Å)')
plt.axvline(x=4226.7, color='purple', linestyle='--', label='Ca I (4226 Å)')
plt.axvline(x=4077.7, color='olive', linestyle='--', label='Sr II (4077 Å)')
plt.axvline(x=4144.0, color='teal', linestyle='--', label='CN (4144 Å)')
plt.axvline(x=4216.0, color='teal', linestyle='--', label='CN (4216 Å)')
plt.axvline(x=4300.0, color='orange', linestyle='--', label='CH (4300 Å)')
plt.axvline(x=4045.8, color='brown', linestyle='--', label='Fe I (4045 Å)')
plt.axvline(x=4271.7, color='brown', linestyle='--', label='Fe I (4271 Å)')
plt.axvline(x=4383.5, color='brown', linestyle='--', label='Fe I (4383 Å)')
plt.title(f"Comparación espectral: estrella 4 vs {tipo_f}")
plt.legend(loc="lower right", ncol=3, fontsize='x-small')
plt.xlim(3900, 5010); plt.ylim(0.0, 1.15)
plt.tight_layout()
plt.show()

#Segunda grafia (lineas de Balmer):
plt.figure(figsize=(12, 6))
plt.plot(w_m_norm, f_m_norm, color='red', alpha=0.6, label=f'MILES ({tipo_f})')
plt.plot(w_prob_norm, f_prob_norm, color='black', label='Estrella 4')

balmer = {
    r'$H\epsilon$ (3970.1 Å)': 3970.1,
    r'$H\delta$ (4101.7 Å)': 4101.7,
    r'$H\gamma$ (4340.5 Å)': 4340.5,
    r'$H\beta$ (4861.3 Å)': 4861.3
}

colors = ['purple', 'blue', 'green', 'darkred']
for i, (name, wave) in enumerate(balmer.items()):
    plt.axvline(x=wave, color=colors[i], linestyle='--', alpha=0.8, label=name)

plt.title(f"Comparación espectral: estrella 4 vs {tipo_f}")
plt.xlabel("Longitud de onda (Å)"); plt.ylabel("Flujo normalizado")
plt.legend(loc="lower right", ncol=3, fontsize='small') 
plt.xlim(3900, 5010); plt.ylim(0.0, 1.15)
plt.tight_layout()
plt.show()