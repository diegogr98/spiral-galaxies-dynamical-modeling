"""
Herramienta de diagnóstico: Muestra valores estadísticos de la extensión SIGMA_SFR.
"""

from pathlib import Path
import numpy as np
from astropy.io import fits

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"


def main():
    fits_file = RAW_DIR / "K0311_gsd6e.fits"
    if not fits_file.exists():
        print(f"Archivo no encontrado: {fits_file}")
        return

    with fits.open(fits_file) as h:
        data = h["SIGMA_SFR"].data

    val = data[np.isfinite(data) & (data > 0)]
    print("Min:", np.min(val), "Max:", np.max(val), "Median:", np.median(val))


if __name__ == "__main__":
    main()
