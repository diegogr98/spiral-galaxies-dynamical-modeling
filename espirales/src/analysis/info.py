"""
Herramienta de diagnóstico: Muestra resumen de estructura HDU en cubos FITS.
"""

from pathlib import Path
from astropy.io import fits

BASE_DIR = Path(__file__).resolve().parent.parent.parent
RAW_DIR = BASE_DIR / "data" / "raw"


def main():
    fits_file = RAW_DIR / "K0311_gsd6e.fits"
    if fits_file.exists():
        with fits.open(fits_file) as hdul:
            hdul.info()
    else:
        print(f"Archivo no encontrado: {fits_file}")


if __name__ == "__main__":
    main()
