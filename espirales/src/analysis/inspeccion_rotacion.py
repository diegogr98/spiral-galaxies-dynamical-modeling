"""
Herramienta de diagnóstico: Diagnóstico de centrado, dispersión y momento angular neto en simulación N-cuerpos.
Entrada: data/raw/disk_000.dat
"""

from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"


def main():
    archivo = RAW_DIR / "disk_000.dat"
    if not archivo.exists():
        print(f"Archivo no encontrado: {archivo}")
        return

    print("--- 1. Primeras 5 líneas crudas ---")
    with open(archivo, "r") as f:
        for _ in range(5):
            print(f.readline().strip())

    print("\n--- 2. Leyendo muestra preliminar (50,000 partículas) ---")
    df_sample = pd.read_csv(
        archivo,
        sep=r"\s+",
        header=None,
        nrows=50000,
        names=["id", "x", "y", "z", "vx", "vy", "vz", "phi_pot"],
    )

    print(f"Total columnas: {list(df_sample.columns)}")

    print("\n--- 3. Diagnóstico de Centrado (Medianas) ---")
    print("Centro espacial (kpc):")
    print(df_sample[["x", "y", "z"]].median().to_dict())

    print("\nVelocidad sistémica residual (km/s):")
    print(df_sample[["vx", "vy", "vz"]].median().to_dict())

    print("\n--- 4. Dispersiones espaciales y cinemáticas (std) ---")
    print("Espesor en ejes (std en kpc):")
    print(df_sample[["x", "y", "z"]].std().to_dict())

    print("\nDispersión de velocidad (std en km/s):")
    print(df_sample[["vx", "vy", "vz"]].std().to_dict())

    # Verificación del eje de rotación: L = r x v
    r = df_sample[["x", "y", "z"]].values
    v = df_sample[["vx", "vy", "vz"]].values
    L = np.cross(r, v)
    L_mean = np.mean(L, axis=0)

    print("\n--- 5. Vector de Momento Angular Neto [Lx, Ly, Lz] (kpc * km/s) ---")
    print(f"Lx: {L_mean[0]:.2e}")
    print(f"Ly: {L_mean[1]:.2e}")
    print(f"Lz: {L_mean[2]:.2e}")


if __name__ == "__main__":
    main()
