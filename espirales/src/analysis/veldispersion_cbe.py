"""
Análisis: Perfiles radiales de dispersión de velocidades estelares para el modelo CBE (zca6e).
Entradas: data/processed/espirales_final_con_bulbo.csv, data/raw/K*_zca6e.fits
Salida: outputs/figures/vel_disp_cbe.png
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent.parent
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from src.utils.fits_tools import calcular_radios, binning_radial

RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
FIGURES_DIR = BASE_DIR / "outputs" / "figures"

# --- CONFIGURACIÓN GENERAL ---
TAMANO_BIN = 2.0  # Tamaño del bin en píxeles
FILTRO_SN = 5  # Umbral de Señal/Ruido
Y_LABEL = r"$\sigma_v \ [\mathrm{km/s}]$"


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    ruta_csv = PROCESSED_DIR / "espirales_final_con_bulbo.csv"
    if not ruta_csv.exists():
        ruta_csv = RAW_DIR / "espirales_final_con_bulbo.csv"

    df_galaxias = pd.read_csv(ruta_csv)
    df_galaxias.columns = df_galaxias.columns.str.strip()

    lista_ids = df_galaxias["califaID"].tolist()
    print(f"Procesando {len(lista_ids)} galaxias para dispersión de velocidades CBE")

    num_galaxias = len(lista_ids)
    columnas = 3
    filas = int(np.ceil(num_galaxias / columnas))

    fig, axes = plt.subplots(filas, columnas, figsize=(12, 3 * filas), sharex=True)
    axes = axes.flatten()

    for i, id_obj in enumerate(lista_ids):
        id_str = f"{id_obj:04d}"
        archivo_fits = RAW_DIR / f"K{id_str}_zca6e.fits"

        fila_galaxia = df_galaxias[df_galaxias["califaID"] == id_obj]
        pa = fila_galaxia["pa"].values[0]
        inc = fila_galaxia["inclination"].values[0]
        re_g = fila_galaxia["radio_bulbo_g"].values[0]
        re_i = fila_galaxia["radio_bulbo_i"].values[0]

        ax = axes[i]

        r_circ, r_el_des, sig_data = calcular_radios(
            archivo_fits, pa, inc, ext_nombre=None, filtro_sn=FILTRO_SN
        )

        if r_circ is not None:
            bc_circ, y_circ = binning_radial(
                r_circ, sig_data, step=TAMANO_BIN, positive_only=True
            )
            bc_el_des, y_el_des = binning_radial(
                r_el_des, sig_data, step=TAMANO_BIN, positive_only=True
            )

            if len(y_circ) > 0:
                ax.plot(
                    bc_circ,
                    y_circ,
                    "o-",
                    color="darkcyan",
                    markersize=3,
                    linewidth=1,
                    label="cbe Circular",
                )
            if len(y_el_des) > 0:
                ax.plot(
                    bc_el_des,
                    y_el_des,
                    "s-",
                    color="indigo",
                    markersize=3,
                    linewidth=1,
                    label="cbe Elíptico",
                )

            # Líneas del radio efectivo del bulbo
            if re_g is not None and not np.isnan(re_g):
                ax.axvline(
                    re_g,
                    color="crimson",
                    linestyle=":",
                    linewidth=1.5,
                    label=f"$R_e$ g ({re_g:.1f}″)",
                )
            if re_i is not None and not np.isnan(re_i):
                ax.axvline(
                    re_i,
                    color="darkorange",
                    linestyle="--",
                    linewidth=1.2,
                    label=f"$R_e$ i ({re_i:.1f}″)",
                )

            ax.set_xlim(0, 40)
            ax.set_ylim(0, 200)
            ax.set_title(f"K{id_str} ($i$={int(inc)}°)", fontsize=10)
            ax.grid(True, alpha=0.2)

            if i >= (filas - 1) * columnas:
                ax.set_xlabel("Radio [píxeles]")
            if i % columnas == 0:
                ax.set_ylabel(Y_LABEL)
            ax.legend(loc="upper right", fontsize=5.5, framealpha=0.6)

        else:
            ax.text(0.5, 0.5, f"FITS K{id_str}\nsin datos", ha="center", va="center")
            ax.set_title(f"K{id_str} (Faltante)")

    for j in range(num_galaxias, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    salida_fig = FIGURES_DIR / "vel_disp_cbe.png"
    plt.savefig(salida_fig, dpi=200)
    plt.close()
    print(f"¡Proceso terminado con éxito! Mosaico guardado como '{salida_fig}'")


if __name__ == "__main__":
    main()
