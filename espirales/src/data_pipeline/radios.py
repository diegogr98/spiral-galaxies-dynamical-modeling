"""
Módulo del pipeline: Extracción de radios efectivos del bulbo (bandas g e i).
Entradas: data/processed/espirales_final.csv (o data/raw/), data/raw/2d_descomposicion_fot.csv
Salida: data/processed/espirales_final_con_bulbo.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"


def main():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    path_espirales = PROCESSED_DIR / "espirales_final.csv"
    if not path_espirales.exists():
        path_espirales = RAW_DIR / "espirales_final.csv"

    print("--- Ejecutando radios.py ---")
    print(f"Cargando muestra de espirales desde: {path_espirales}")
    df_espirales = pd.read_csv(path_espirales)
    df_espirales.columns = df_espirales.columns.str.strip()

    path_photdec = RAW_DIR / "2d_descomposicion_fot.csv"
    print(f"Cargando descomposición fotométrica 2D desde: {path_photdec}")
    df_photdec = pd.read_csv(
        path_photdec, comment="#", header=None, skipinitialspace=True
    )

    # Columnas: 0 -> CALIFAID, 14 -> Re_g, 126 -> Re_i
    df_radios = df_photdec[[0, 14, 126]].copy()
    df_radios.columns = ["CALIFAID", "radio_bulbo_g", "radio_bulbo_i"]

    df_espirales["califaID"] = df_espirales["califaID"].astype(int)
    df_radios["CALIFAID"] = df_radios["CALIFAID"].astype(int)

    print("Efectuando la unión de los catálogos...")
    df_final = pd.merge(
        df_espirales,
        df_radios,
        left_on="califaID",
        right_on="CALIFAID",
        how="left",
    )
    df_final = df_final.drop(columns=["CALIFAID"])

    # Limpieza de valores nulos astronómicos (-999.0)
    df_final["radio_bulbo_g"] = df_final["radio_bulbo_g"].replace(-999.0, np.nan)
    df_final["radio_bulbo_i"] = df_final["radio_bulbo_i"].replace(-999.0, np.nan)

    salida = PROCESSED_DIR / "espirales_final_con_bulbo.csv"
    df_final.to_csv(salida, index=False)

    print(f"\n¡Fusión completada exitosamente!")
    print(f"Archivo guardado como: {salida}")
    print(f"Columnas resultantes: {list(df_final.columns)}")
    print("\nMuestra de los primeros datos combinados:")
    print(df_final[["califaID", "radio_bulbo_g", "radio_bulbo_i"]].dropna().head())


if __name__ == "__main__":
    main()
