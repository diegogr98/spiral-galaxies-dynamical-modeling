"""
Módulo de utilidades para lectura, análisis y extracción de perfiles en archivos FITS de CALIFA.
"""

from pathlib import Path
import numpy as np
from astropy.io import fits


def calcular_radios(archivo, pa_deg, inc_deg, ext_nombre=None, filtro_sn=5):
    """
    Calcula mallas radiales (circular y elíptica desinclinada) y extrae datos de un archivo FITS.

    Parámetros:
        archivo (str or Path): Ruta al archivo FITS.
        pa_deg (float): Ángulo de posición en grados.
        inc_deg (float): Inclinación galáctica en grados.
        ext_nombre (str, opcional): Nombre de la extensión FITS (ej: 'L_5635', 'V_D', 'SIG_V').
                                    Si es None, busca 'V_D' o 'SIG_V'.
        filtro_sn (float, opcional): Umbral mínimo de Señal/Ruido (S/N). Por defecto 5.

    Retorna:
        tuple: (r_circular, r_elip_desinc, data_map) aplanados en 1D, o (None, None, None) si falla.
    """
    archivo = Path(archivo)
    if not archivo.exists():
        return None, None, None

    phi = np.radians(pa_deg + 90.0)
    cos_i = np.cos(np.radians(inc_deg))

    with fits.open(archivo) as hdul:
        if ext_nombre is None:
            ext = "V_D" if "V_D" in hdul else "SIG_V"
        else:
            ext = ext_nombre

        if ext not in hdul:
            print(f"Advertencia: No se encontró la extensión {ext} en {archivo.name}")
            return None, None, None

        data_map = hdul[ext].data.copy()
        sn = hdul["SN"].data if "SN" in hdul else None
        header = hdul[ext].header

        # Máscara usando el canal de Señal/Ruido
        if sn is not None and filtro_sn is not None:
            data_map[sn < filtro_sn] = np.nan

        # Centro astronómico (CRPIX) tomado del header del FITS
        xc = header.get("CRPIX1", data_map.shape[1] / 2.0)
        yc = header.get("CRPIX2", data_map.shape[0] / 2.0)

        # Crear mallas de coordenadas espaciales
        y, x = np.indices(data_map.shape)
        xt, yt = x - xc, y - yc

        # Modelado de geometrías radiales (Circular vs Elíptica Desinclinada)
        r_circular = np.sqrt(xt**2 + yt**2)
        x_prime = xt * np.cos(phi) + yt * np.sin(phi)
        y_prime = -xt * np.sin(phi) + yt * np.cos(phi)
        r_elip_desinc = np.sqrt(x_prime**2 + (y_prime / cos_i)**2)

        return r_circular.flatten(), r_elip_desinc.flatten(), data_map.flatten()


def binning_radial(r, data, step, positive_only=False):
    """
    Realiza un perfil radial azimutal calculando la mediana en anillos de ancho `step`.

    Parámetros:
        r (ndarray): Array de radios (píxeles o arcsec).
        data (ndarray): Array de valores a promediar/medir.
        step (float): Ancho del bin radial.
        positive_only (bool): Si True, solo incluye bins con mediana > 0 (útil para dispersión).

    Retorna:
        tuple: (bin_centers, bin_values)
    """
    if len(r) == 0 or np.all(np.isnan(r)):
        return np.array([]), np.array([])

    max_r = np.nanmax(r)
    if np.isnan(max_r) or max_r <= 0:
        return np.array([]), np.array([])

    bins = np.arange(0, max_r, step)
    bc, bv = [], []
    for i in range(len(bins) - 1):
        mask = (r >= bins[i]) & (r < bins[i + 1])
        if np.any(mask) and not np.all(np.isnan(data[mask])):
            mediana = np.nanmedian(data[mask])
            if np.isnan(mediana):
                continue
            if positive_only and mediana <= 0:
                continue
            bc.append((bins[i] + bins[i + 1]) / 2.0)
            bv.append(mediana)

    return np.array(bc)[:len(bv)], np.array(bv)
