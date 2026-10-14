"""
Módulo del pipeline: Filtro cinemático y morfológico para selección de galaxias espirales puras.
Entradas: data/raw/gme.fits, data/raw/cbe.fits, data/raw/growth_curves.csv, data/raw/morphology.csv
Salida: data/processed/espirales_final.csv
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from astropy.table import Table

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"


def main():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    print("--- Ejecutando filtro_incl.py ---")
    print(f"Leyendo catálogos crudos desde: {RAW_DIR}")

    # 1. Carga de catálogos crudos
    df_gme = Table.read(RAW_DIR / "gme.fits").to_pandas()
    df_cbe = Table.read(RAW_DIR / "cbe.fits").to_pandas()
    df_geo = pd.read_csv(RAW_DIR / "growth_curves.csv", comment="#", header=None).rename(
        columns={0: "califaID", 1: "pa", 2: "ba"}
    )
    df_morf = pd.read_csv(RAW_DIR / "morphology.csv", comment="#", header=None).rename(
        columns={0: "califaID", 5: "hubtyp", 6: "subtyp"}
    )

    # Limpieza y unión
    df_gme["califaID"] = df_gme["califaID"].astype(int)
    df_cbe["califaID"] = df_cbe["califaID"].astype(int)

    df_pyc = pd.merge(df_gme, df_cbe, on="califaID", suffixes=("_gme", "_cbe"))
    df_total = pd.merge(df_pyc, df_morf[["califaID", "hubtyp", "subtyp"]], on="califaID")
    df_total = pd.merge(df_total, df_geo[["califaID", "pa", "ba"]], on="califaID")

    # Calculando inclinaciones
    q0 = 0.0
    cos_i2 = (df_total["ba"] ** 2 - q0**2) / (1.0 - q0**2)
    df_total["inclination"] = np.degrees(np.arccos(np.sqrt(cos_i2.clip(lower=0))))

    # Máscara estricta: i < 30, Sigma > 80
    mask_sigma = df_total["sigma_V_gme"] > 80
    mask_inc = df_total["inclination"] < 30
    mask_espiral = df_total["hubtyp"].astype(str).str.contains("S", na=False)
    mask_sin_s0 = ~df_total["subtyp"].astype(str).str.contains("0", na=False)

    df_espirales = df_total[mask_sigma & mask_inc & mask_espiral & mask_sin_s0].copy()

    # Guardar CSV procesado
    columnas_finales = ["califaID", "hubtyp", "subtyp", "inclination", "pa", "sigma_V_gme"]
    salida = PROCESSED_DIR / "espirales_final.csv"
    df_espirales[columnas_finales].to_csv(salida, index=False)

    print(f"Éxito: Se encontraron {len(df_espirales)} espirales puras.")
    print(f"Archivo guardado como: {salida}")


if __name__ == "__main__":
    main()
