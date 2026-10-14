"""
Análisis: Exploración del espacio de parámetros de orientación (PA e inclinación) y degeneración de multiplicidad.
Entrada: data/raw/disk_000.dat
Salida: outputs/figures/rotacion_multiplicidad.png
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from src.utils.geometria import rotar_proyeccion

RAW_DIR = BASE_DIR / "data" / "raw"
FIGURES_DIR = BASE_DIR / "outputs" / "figures"


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Cargar archivo completo
    archivo = RAW_DIR / "disk_000.dat"
    print(f"Cargando archivo de partículas desde: {archivo}")

    df = pd.read_csv(
        archivo,
        sep=r"\s+",
        header=None,
        names=["id", "x", "y", "z", "vx", "vy", "vz", "phi_pot"],
    )
    print(f"Total partículas leídas: {len(df):,}")

    # 2. Centrado estricto (restar medianas)
    pos0 = np.vstack([df["x"].values, df["y"].values, df["z"].values])
    vel0 = np.vstack([df["vx"].values, df["vy"].values, df["vz"].values])

    pos0 -= np.median(pos0, axis=1, keepdims=True)
    vel0 -= np.median(vel0, axis=1, keepdims=True)

    # 3. Parámetros de prueba (CALIFA 813)
    pa_base = 98.53
    inc_base = 19.42

    # Combinaciones de multiplicidad para comparar
    casos = [
        (pa_base, inc_base, f"PA={pa_base:.1f}°, inc=+{inc_base:.1f}°"),
        (
            pa_base + 180.0,
            inc_base,
            f"PA={pa_base+180:.1f}°, inc=+{inc_base:.1f}° (PA+180)",
        ),
        (pa_base, -inc_base, f"PA={pa_base:.1f}°, inc=-{inc_base:.1f}° (-inc)"),
        (
            pa_base + 180.0,
            -inc_base,
            f"PA={pa_base+180:.1f}°, inc=-{inc_base:.1f}° (Ambos)",
        ),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    rango = 25.0  # kpc

    for ax, (pa_val, inc_val, titulo) in zip(axes, casos):
        p_rot, v_rot = rotar_proyeccion(pos0, vel0, pa_val, inc_val)
        x_proj, y_proj = p_rot[0], p_rot[1]
        v_los = v_rot[2]

        # Mapa 2D de velocidad en línea de visión (v_los)
        hb = ax.hexbin(
            x_proj,
            y_proj,
            C=v_los,
            reduce_C_function=np.median,
            gridsize=90,
            cmap="coolwarm",
            extent=[-rango, rango, -rango, rango],
            vmin=-200,
            vmax=200,
        )
        ax.set_title(titulo, fontsize=11)
        ax.set_xlabel("X proyectado (kpc)")
        ax.set_ylabel("Y proyectado (kpc)")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)
        fig.colorbar(hb, ax=ax, label=r"$V_{\rm los}$ (km/s)")

    plt.tight_layout()
    salida_fig = FIGURES_DIR / "rotacion_multiplicidad.png"
    plt.savefig(salida_fig, dpi=200)
    plt.close()
    print(f"Gráfica guardada en {salida_fig}")


if __name__ == "__main__":
    main()