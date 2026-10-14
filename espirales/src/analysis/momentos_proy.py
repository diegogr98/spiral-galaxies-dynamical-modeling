"""
Análisis: Generación de momentos cinemáticos y mapas sintéticos de spaxels a partir de simulaciones N-cuerpos.
Entrada: data/raw/disk_000.dat
Salida: outputs/figures/mapas_sinteticos_spaxels.png
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d

BASE_DIR = Path(__file__).resolve().parent.parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from src.utils.geometria import rotar_proyeccion

RAW_DIR = BASE_DIR / "data" / "raw"
FIGURES_DIR = BASE_DIR / "outputs" / "figures"


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # --- 1. Cargar archivo de partículas ---
    archivo = RAW_DIR / "disk_000.dat"
    print(f"Cargando archivo de partículas desde: {archivo}")
    df = pd.read_csv(
        archivo,
        sep=r"\s+",
        header=None,
        names=["id", "x", "y", "z", "vx", "vy", "vz", "phi_pot"],
    )

    pos0 = np.vstack([df["x"].values, df["y"].values, df["z"].values])
    vel0 = np.vstack([df["vx"].values, df["vy"].values, df["vz"].values])

    pos0 -= np.median(pos0, axis=1, keepdims=True)
    vel0 -= np.median(vel0, axis=1, keepdims=True)

    # --- 2. Parámetros geométricos y escala de CALIFA (CALIFA 813) ---
    pa_val = 98.53
    inc_val = 19.42

    distancia_mpc = 65.0
    kpc_por_arcsec = (distancia_mpc * 1e3) / 206265.0

    # --- 3. Rotación y proyección ---
    p_rot, v_rot = rotar_proyeccion(pos0, vel0, pa_val, inc_val)

    # Reescalar posiciones a segundos de arco
    x_arcsec = p_rot[0] / kpc_por_arcsec
    y_arcsec = p_rot[1] / kpc_por_arcsec
    v_los = v_rot[2]

    # --- 4. Grilla de spaxels sintéticos con FOV ampliado a 70.0'' ---
    fov_arcsec = 70.0
    spaxel_size = 1.0
    bins = np.arange(-fov_arcsec, fov_arcsec + spaxel_size, spaxel_size)

    # Binned statistic: Velocidad media y Dispersión por celda
    v_mean, x_edges, y_edges, _ = binned_statistic_2d(
        x_arcsec, y_arcsec, v_los, statistic="mean", bins=[bins, bins]
    )

    sigma_map, _, _, bin_counts = binned_statistic_2d(
        x_arcsec, y_arcsec, v_los, statistic="std", bins=[bins, bins]
    )

    # Corrección de Bug Crítico: pasar v_los en lugar de None
    conteo, _, _, _ = binned_statistic_2d(
        x_arcsec, y_arcsec, v_los, statistic="count", bins=[bins, bins]
    )

    # Enmascarar celdas con conteo < 15
    mask = conteo < 15
    v_mean[mask] = np.nan
    sigma_map[mask] = np.nan

    # --- 5. Visualización de los mapas sintéticos ---
    extent = [-fov_arcsec, fov_arcsec, -fov_arcsec, fov_arcsec]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    # Mapa de Velocidad Media (V_los)
    im1 = ax1.imshow(
        v_mean.T, origin="lower", extent=extent, cmap="coolwarm", vmin=-150, vmax=150
    )
    ax1.set_title(r"Mapa Sintético $V_{\rm los}$ (Media)", fontsize=13)
    ax1.set_xlabel(r"$\Delta \alpha$ (arcsec)")
    ax1.set_ylabel(r"$\Delta \delta$ (arcsec)")
    fig.colorbar(im1, ax=ax1, label=r"$V_{\rm los}$ (km/s)", fraction=0.046, pad=0.04)

    # Mapa de Dispersión de Velocidades (sigma_los)
    im2 = ax2.imshow(
        sigma_map.T, origin="lower", extent=extent, cmap="viridis", vmin=0, vmax=120
    )
    ax2.set_title(r"Mapa Sintético $\sigma_{\rm los}$ (Dispersión)", fontsize=13)
    ax2.set_xlabel(r"$\Delta \alpha$ (arcsec)")
    fig.colorbar(im2, ax=ax2, label=r"$\sigma_{\rm los}$ (km/s)", fraction=0.046, pad=0.04)

    for ax in (ax1, ax2):
        ax.axhline(0, color="grey", linestyle="--", alpha=0.5)
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)

    plt.tight_layout()
    salida_fig = FIGURES_DIR / "mapas_sinteticos_spaxels.png"
    plt.savefig(salida_fig, dpi=200)
    plt.close()
    print(f"Mapas sintéticos generados en {salida_fig}")


if __name__ == "__main__":
    main()