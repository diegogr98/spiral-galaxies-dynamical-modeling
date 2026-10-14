"""
Módulo del pipeline: Consolidación de parámetros geométricos, descomposiciones 2D y masas estelares.
Entradas: data/processed/espirales_final_con_bulbo.csv, data/raw/2d_descomposicion_fot.csv, data/raw/componentes_Gal.csv
Salida: data/processed/parametros.csv
"""

from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"


def main():
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    path_lista = PROCESSED_DIR / "espirales_final_con_bulbo.csv"
    if not path_lista.exists():
        path_lista = PROCESSED_DIR / "espirales_final.csv"
        if not path_lista.exists():
            path_lista = RAW_DIR / "espirales_final.csv"

    path_descomp2d = RAW_DIR / "2d_descomposicion_fot.csv"
    path_masas = RAW_DIR / "componentes_Gal.csv"

    print("--- Ejecutando parametros.py ---")
    print(f"Cargando lista base desde: {path_lista}")
    # 1. Cargar df_lista
    df_lista = pd.read_csv(path_lista)[["califaID", "inclination", "pa"]]
    df_lista["califaID"] = df_lista["califaID"].astype(int)

    # 2. Cargar df_descomp2d y reemplazar los -999.0 por NaN
    print(f"Cargando descomposición fotométrica desde: {path_descomp2d}")
    df_descomp2d = pd.read_csv(
        path_descomp2d,
        skiprows=188,
        header=None,
        usecols=[0, 70, 81, 85, 87, 89],
        names=["califaID", "Re", "hi_r", "PAD_r", "ho_r", "r_break"],
    )
    df_descomp2d["califaID"] = df_descomp2d["califaID"].astype(int)
    df_descomp2d.replace(-999.0, np.nan, inplace=True)

    # 3. Cargar df_masas
    print(f"Cargando masas estelares por componente desde: {path_masas}")
    df_masas = pd.read_csv(
        path_masas,
        skiprows=54,
        header=None,
        usecols=[0, 2, 4, 6],
        names=["califaID", "M_bulgue", "M_bar", "M_disk"],
    )
    df_masas["califaID"] = df_masas["califaID"].astype(int)
    df_masas.replace(-999.0, np.nan, inplace=True)

    # 4. Unir catálogos
    df_final = df_lista.merge(df_descomp2d, on="califaID", how="left").merge(
        df_masas, on="califaID", how="left"
    )

    salida = PROCESSED_DIR / "parametros.csv"
    df_final.to_csv(salida, index=False)

    print("Catálogo Final Unido ------")
    print(df_final)
    print(f"\nArchivo guardado como: {salida}")


if __name__ == "__main__":
    main()