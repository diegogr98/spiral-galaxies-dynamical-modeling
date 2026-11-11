"""
Análisis: Generación de momentos cinemáticos y mapas sintéticos de spaxels a partir de simulaciones N-cuerpos.
Entrada: disco_684.dat (N-cuerpos), data/raw/gme.fits, data/processed/parametros_con_gas.csv
Salida: outputs/figures/mapas_sinteticos_spaxels_684.png
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from astropy.table import Table

BASE_DIR = Path(__file__).resolve().parent.parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from src.utils.geometria import rotar_proyeccion

RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
FIGURES_DIR = BASE_DIR / "outputs" / "figures"


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Definir la galaxia objetivo a la 684
    target_id = 684

    # --- 1. Extraer Parámetros Geométricos (CSV) y Distancia (FITS) ---
    print(f"Extrayendo parámetros observacionales para CALIFA {target_id}...")
    
    # Extraer PA e Inclinación
    csv_path = PROCESSED_DIR / "parametros_con_gas.csv"
    df_params = pd.read_csv(csv_path)
    galaxia_csv = df_params[df_params["califaID"] == target_id].iloc[0]

    pa_val = galaxia_csv["pa_calc"] + 180.0
    inc_val = galaxia_csv["inclination"]
    
    # Extraer Distancia
    fits_path = RAW_DIR / "gme.fits"
    df_gme = Table.read(fits_path).to_pandas()
    galaxia_fits = df_gme[df_gme["califaID"] == target_id].iloc[0]
    distancia_mpc = galaxia_fits["DL_Mpc"]
    d_kpc = distancia_mpc * 1000.0  # Distancia al observador en kpc
    
    print(f" > PA: {pa_val:.2f}°, Inclinación: {inc_val:.2f}°, Distancia: {distancia_mpc:.2f} Mpc")

    # --- 2. Cargar archivo de partículas ---
    # Usando la ruta exacta de la nueva simulación
    archivo = Path("/home/diego/Documents/tesis/Galactic_dynamics_Thesis/espirales/simulacion_galaxia/AM2322A_MakeDisk.isolated/disco_684.dat")
    print(f"\nCargando archivo de partículas desde: {archivo}")
    
    # El archivo .dat generado por el C solo tiene 6 columnas
    df = pd.read_csv(
        archivo,
        sep=r"\s+",
        header=None,
        names=["x", "y", "z", "vx", "vy", "vz"],
    )

    pos0 = np.vstack([df["x"].values, df["y"].values, df["z"].values])
    vel0 = np.vstack([df["vx"].values, df["vy"].values, df["vz"].values])

    pos0 -= np.median(pos0, axis=1, keepdims=True)
    vel0 -= np.median(vel0, axis=1, keepdims=True)

    # --- 3. Rotación y proyección con corrección de perspectiva ---
    kpc_por_arcsec = d_kpc / 206265.0
    
    # Rotación 3D
    p_rot, v_rot = rotar_proyeccion(pos0, vel0, pa_val, inc_val)
    x_proj, y_proj, z_proj = p_rot[0], p_rot[1], p_rot[2]
    vx, vy, vz = v_rot[0], v_rot[1], v_rot[2]

    # Calcular el vector unitario de la línea de visión (r_hat) para cada partícula
    rx_los = x_proj
    ry_los = y_proj
    rz_los = z_proj - d_kpc 
    
    norm_r = np.sqrt(rx_los**2 + ry_los**2 + rz_los**2)
    r_hat_x = rx_los / norm_r
    r_hat_y = ry_los / norm_r
    r_hat_z = rz_los / norm_r
    
    # Producto punto para obtener la velocidad en la línea de visión exacta
    v_los = (vx * r_hat_x) + (vy * r_hat_y) + (vz * r_hat_z)

    # Reescalar posiciones a segundos de arco para la grilla
    x_arcsec = x_proj / kpc_por_arcsec
    y_arcsec = y_proj / kpc_por_arcsec

    # --- 4. Grilla de spaxels sintéticos con FOV ampliado a 45.0'' ---
    fov_arcsec = 45.0
    spaxel_size = 2.0
    bins = np.arange(-fov_arcsec, fov_arcsec + spaxel_size, spaxel_size)

    # Binned statistic: Velocidad media y Dispersión por celda
    v_mean, x_edges, y_edges, _ = binned_statistic_2d(
        x_arcsec, y_arcsec, v_los, statistic="mean", bins=[bins, bins]
    )

    sigma_map, _, _, bin_counts = binned_statistic_2d(
        x_arcsec, y_arcsec, v_los, statistic="std", bins=[bins, bins]
    )

    conteo, _, _, _ = binned_statistic_2d(
        x_arcsec, y_arcsec, v_los, statistic="count", bins=[bins, bins]
    )

    # Enmascarar celdas con conteo < 8
    mask = conteo < 8
    v_mean[mask] = np.nan
    sigma_map[mask] = np.nan

    # --- 5. Visualización de los mapas sintéticos ---
    extent = [-fov_arcsec, fov_arcsec, -fov_arcsec, fov_arcsec]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    # Mapa de Velocidad Media (V_los)
    im1 = ax1.imshow(
        v_mean.T, origin="lower", extent=extent, cmap="RdBu_r", vmin=-200, vmax=200
    )
    ax1.set_title(rf"Mapa Sintético $V_{{\rm los}}$ (Media) - CALIFA {target_id}", fontsize=13)
    ax1.set_xlabel(r"$\Delta \alpha$ (arcsec)")
    ax1.set_ylabel(r"$\Delta \delta$ (arcsec)")
    fig.colorbar(im1, ax=ax1, label=r"$V_{\rm los}$ (km/s)", fraction=0.046, pad=0.04)

    # Mapa de Dispersión de Velocidades (sigma_los)
    im2 = ax2.imshow(
        sigma_map.T, origin="lower", extent=extent, cmap="inferno_r", vmin=0, vmax=220
    )
    ax2.set_title(rf"Mapa Sintético $\sigma_{{\rm los}}$ (Dispersión) - CALIFA {target_id}", fontsize=13)
    ax2.set_xlabel(r"$\Delta \alpha$ (arcsec)")
    fig.colorbar(im2, ax=ax2, label=r"$\sigma_{\rm los}$ (km/s)", fraction=0.046, pad=0.04)

    for ax in (ax1, ax2):
        ax.axhline(0, color="grey", linestyle="--", alpha=0.5)
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)

    plt.tight_layout()
    salida_fig = FIGURES_DIR / f"mapas_sinteticos_spaxels_{target_id}.png"
    plt.savefig(salida_fig, dpi=200)
    plt.close()
    print(f"\nMapas sintéticos generados en {salida_fig}")


if __name__ == "__main__":
    main()