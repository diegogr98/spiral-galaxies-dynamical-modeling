"""
Herramienta de diagnóstico: Inspección preliminar de columnas y filas en gme.fits y cbe.fits.
"""

from pathlib import Path
import pandas as pd
from astropy.table import Table

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"


def main():
    try:
        df_gme = Table.read(RAW_DIR / "gme.fits").to_pandas()
        df_cbe = Table.read(RAW_DIR / "cbe.fits").to_pandas()

        print("--- COLUMNAS ENCONTRADAS EN GMe ---")
        print(df_gme.columns.tolist())

        print("\n--- PRIMERAS FILAS PARA VER VALORES ---")
        print(df_gme.head(20))

    except Exception as e:
        print(f"Error al leer los archivos: {e}")


if __name__ == "__main__":
    main()
