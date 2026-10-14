"""
Módulo del pipeline: Cálculo de masa y escala de gas molecular/atómico mediante inversión de Kennicutt-Schmidt.
Entradas: data/processed/parametros.csv, data/raw/K*_gsd6e.fits
Salidas: data/processed/parametros_con_gas.csv, outputs/figures/perfiles_gas_final.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
FIGURES_DIR = BASE_DIR / "outputs" / "figures"

# 1. Constantes empíricas de Schmidt-Kennicutt (Kennicutt 1998)
A_ks = 2.5e-4
N_ks = 1.4

# 2. Rangos de ajuste (arcsec) donde el disco presenta decaimiento exponencial continuo
FIT_RANGES = {
    311: (2.0, 9.0),
    684: (19.0, 33.0),
    740: (5.0, 25.0),
    813: (16.0, 33.0),
    872: (4.0, 12.0),
    915: (7.0, 25.0),
}


def main():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    csv_path = PROCESSED_DIR / "parametros.csv"
    if not csv_path.exists():
        csv_path = RAW_DIR / "parametros.csv"

    print("--- Ejecutando calcular_gas.py ---")
    print(f"Cargando parámetros desde: {csv_path}")
    df = pd.read_csv(csv_path)
    df["pa_calc"] = df["PAD_r"].fillna(df["pa"])

    results = []
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, (idx, row) in enumerate(df.iterrows()):
        cid = int(row["califaID"])
        fits_file = RAW_DIR / f"K{cid:04d}_gsd6e.fits"

        if not fits_file.exists():
            print(f"Archivo no encontrado: {fits_file}")
            results.append({
                "califaID": cid,
                "DL_Mpc": np.nan,
                "h_gas_arcsec": np.nan,
                "Sigma0_gas": np.nan,
                "M_gas_1e10": np.nan,
                "M_gas_Msun": np.nan,
            })
            continue

        with fits.open(fits_file) as hdul:
            hdr = hdul["SIGMA_SFR"].header
            sfr_pycasso = hdul["SIGMA_SFR"].data.astype(float)
            dl_mpc = hdul[0].header.get("DL_MPC", None)

        x0 = hdr["CRPIX1"] - 1.0
        y0 = hdr["CRPIX2"] - 1.0
        inc = np.radians(row["inclination"])
        pa = np.radians(row["pa_calc"] - 90.0)

        # Conversión formal de unidades y desproyección intrínseca
        # M_sun / Gyr / pc^2 -> M_sun / yr / kpc^2 corregido a visión frontal (face-on)
        valid_mask = np.isfinite(sfr_pycasso) & (sfr_pycasso > 0)
        sfr_kpc2_faceon = (sfr_pycasso * 1e-3) * np.cos(inc)

        # Inversión de Schmidt-Kennicutt en el plano físico del disco -> Sigma_gas real [M_sun / pc^2]
        sigma_gas = np.full_like(sfr_pycasso, np.nan)
        sigma_gas[valid_mask] = (sfr_kpc2_faceon[valid_mask] / A_ks) ** (1.0 / N_ks)

        # Transformación a coordenadas intrínsecas desproyectadas en el disco
        ny, nx = sigma_gas.shape
        y, x = np.mgrid[:ny, :nx]
        dx = x - x0
        dy = y - y0
        x_rot = dx * np.cos(pa) + dy * np.sin(pa)
        y_rot = -dx * np.sin(pa) + dy * np.cos(pa)
        r_gal = np.sqrt(x_rot**2 + (y_rot / np.cos(inc))**2)

        # Perfil radial azimutal promedio (anillos de 1'')
        r_bins = np.arange(0, 36.0, 1.0)
        r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
        profile_gas = []

        for r_in, r_out in zip(r_bins[:-1], r_bins[1:]):
            ring = (r_gal >= r_in) & (r_gal < r_out) & np.isfinite(sigma_gas)
            if np.sum(ring) >= 5:
                profile_gas.append(np.nanmean(sigma_gas[ring]))
            else:
                profile_gas.append(np.nan)
        profile_gas = np.array(profile_gas)

        # Ajuste del disco exponencial sobre el rango asignado
        r_min_fit, r_max_fit = FIT_RANGES.get(cid, (5.0, 25.0))
        fit_mask = (
            (r_centers >= r_min_fit)
            & (r_centers <= r_max_fit)
            & np.isfinite(profile_gas)
            & (profile_gas > 0)
        )

        if np.sum(fit_mask) >= 3:
            m, b = np.polyfit(r_centers[fit_mask], np.log(profile_gas[fit_mask]), deg=1)
            h_gas = -1.0 / m if m < 0 else np.nan
            sigma0 = np.exp(b) if m < 0 else np.nan
        else:
            h_gas, sigma0 = np.nan, np.nan

        # Integración de la masa total de gas utilizando DL_MPC oficial
        if dl_mpc is not None and dl_mpc > 0:
            pc_per_arcsec = (dl_mpc * 1e6) / 206265.0
            # Área real desproyectada del spaxel en el plano físico de la galaxia
            area_spaxel_pc2 = (pc_per_arcsec**2) / np.cos(inc)
            m_gas_total = np.nansum(sigma_gas) * area_spaxel_pc2
            m_gas_1e10 = m_gas_total / 1e10
        else:
            m_gas_total = np.nan
            m_gas_1e10 = np.nan

        results.append({
            "califaID": cid,
            "DL_Mpc": dl_mpc,
            "h_gas_arcsec": h_gas,
            "Sigma0_gas": sigma0,
            "M_gas_1e10": m_gas_1e10,
            "M_gas_Msun": m_gas_total,
        })

        # Gráfica diagnóstica
        if i < len(axes):
            ax = axes[i]
            ax.plot(r_centers, profile_gas, "ko-", markersize=4, label=r"$\Sigma_{\rm gas}$")
            if np.isfinite(h_gas):
                r_fit = np.linspace(r_min_fit, r_max_fit, 50)
                ax.plot(
                    r_fit,
                    sigma0 * np.exp(-r_fit / h_gas),
                    "r--",
                    lw=2,
                    label=f'h={h_gas:.1f}"',
                )
            ax.set_yscale("log")
            ax.set_title(
                f"CALIFA {cid} ($M_{{\\rm gas}}={m_gas_1e10:.3f}\\times 10^{{10}} M_\\odot$)"
            )
            ax.set_xlabel("Radio (arcsec)")
            ax.set_ylabel(r"$\Sigma_{\rm gas}$ [$M_\odot/{\rm pc}^2$]")
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=8)

        print(
            f"CALIFA {cid:03d} | h_gas = {h_gas:5.2f}'' | M_gas = {m_gas_1e10:.3f} x 10^10 M_sun ({m_gas_total:.2e} M_sun)"
        )

    plt.tight_layout()
    salida_fig = FIGURES_DIR / "perfiles_gas_final.png"
    plt.savefig(salida_fig, dpi=200)
    plt.close()
    print(f"\nGráfica guardada en {salida_fig}")

    # Guardar catálogo actualizado
    df_gas = pd.DataFrame(results)
    cols_to_add = ["califaID", "h_gas_arcsec", "Sigma0_gas", "M_gas_1e10"]
    # Si las columnas ya existen en df, eliminarlas antes del merge para evitar sufijos _x, _y
    cols_existentes = [c for c in cols_to_add[1:] if c in df.columns]
    if cols_existentes:
        df = df.drop(columns=cols_existentes)

    df_final = df.merge(df_gas[cols_to_add], on="califaID", how="left")
    salida_csv = PROCESSED_DIR / "parametros_con_gas.csv"
    df_final.to_csv(salida_csv, index=False)
    print(f"Archivo final guardado en {salida_csv}")


if __name__ == "__main__":
    main()