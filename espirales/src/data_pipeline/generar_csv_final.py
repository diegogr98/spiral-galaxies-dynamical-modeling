"""
Módulo del pipeline: Consolidación de parámetros finales para la simulación.
Une datos tabulares procesados (incluyendo gas) con el catálogo FITS (gme.fits) de CALIFA
y calcula la escala de longitud efectiva (h_optimo) para galaxias de dos discos.

Entradas:
    - data/processed/parametros_con_gas.csv
    - data/raw/gme.fits
Salidas:
    - data/processed/parametros_finales_simulacion.csv
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table

# Paso 1: Configuración de Rutas
BASE_DIR = Path(__file__).resolve().parent.parent
if not (BASE_DIR / "data").exists() and (BASE_DIR.parent / "data").exists():
    BASE_DIR = BASE_DIR.parent

sys.path.insert(0, str(BASE_DIR))

from src.data_pipeline.calcular_h_optimo import calculate_effective_h


def calcular_h_fila(row: pd.Series) -> float:
    """
    Calcula h_optimo para una galaxia si posee parámetros de dos discos (hi_r, ho_r, r_break).
    En caso contrario, retorna np.nan.
    """
    hi = row.get("hi_r")
    ho = row.get("ho_r")
    r_break = row.get("r_break")

    if pd.notna(hi) and pd.notna(ho) and pd.notna(r_break):
        rmax_por_galaxia = {684: 45.0, 813: 50.0}
        cid = int(row["califaID"]) if pd.notna(row.get("califaID")) else -1
        r_max = rmax_por_galaxia.get(cid, 40.0)

        resultado = calculate_effective_h(
            float(hi), float(ho), float(r_break), r_max=r_max
        )
        return float(resultado[0]) if isinstance(resultado, tuple) else float(resultado)

    return np.nan


def main():
    print("--- Iniciando consolidación de parámetros finales (generar_csv_final.py) ---")

    ruta_csv_entrada = BASE_DIR / "data" / "processed" / "parametros_con_gas.csv"
    ruta_fits = BASE_DIR / "data" / "raw" / "gme.fits"
    ruta_csv_salida = BASE_DIR / "data" / "processed" / "parametros_finales_simulacion.csv"

    # Paso 1: Carga del CSV principal y limpieza manual
    print(f"[Paso 1/5] Cargando archivo CSV principal: {ruta_csv_entrada}")
    try:
        if not ruta_csv_entrada.exists():
            raise FileNotFoundError(f"No se encontró el archivo CSV: {ruta_csv_entrada}")
        df_principal = pd.read_csv(ruta_csv_entrada)
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Ocurrió un problema al leer '{ruta_csv_entrada}': {e}")
        sys.exit(1)

    # Limpieza de espacios en blanco en columnas numéricas si existen
    df_principal = df_principal.replace(r"^\s*$", np.nan, regex=True)

    galaxias_excluidas = [802, 925]
    filas_antes = len(df_principal)
    df_principal = df_principal[~df_principal["califaID"].isin(galaxias_excluidas)].copy()
    print(
        f"  -> Galaxias excluidas (califaID {galaxias_excluidas}). "
        f"Filas restantes: {len(df_principal)} (de {filas_antes})."
    )

    # Paso 2: Extracción de Masa Estelar desde catálogo FITS
    print(f"[Paso 2/5] Leyendo catálogo FITS para extraer M_star: {ruta_fits}")
    try:
        if not ruta_fits.exists():
            raise FileNotFoundError(f"No se encontró el catálogo FITS: {ruta_fits}")
        df_gme = Table.read(ruta_fits).to_pandas()
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Ocurrió un problema al leer el archivo FITS '{ruta_fits}': {e}")
        sys.exit(1)

    df_mstar = df_gme[["califaID", "M_star"]].drop_duplicates(subset=["califaID"])
    df_mstar["M_star"] = df_mstar["M_star"] / 1e10

    df_consolidado = pd.merge(df_principal, df_mstar, on="califaID", how="left")
    print("  -> Cruce (left merge) de M_star completado exitosamente.")

    # Paso 3: Cálculo de h_optimo para galaxias de dos discos
    print("[Paso 3/5] Calculando h_optimo para galaxias con doble disco...")
    df_consolidado["h_optimo"] = df_consolidado.apply(calcular_h_fila, axis=1)
    calculados = df_consolidado["h_optimo"].notna().sum()
    print(f"  -> h_optimo calculado para {calculados} galaxia(s); el resto asignado a NaN.")

    # Paso 4: Filtrado y ordenamiento estricto de columnas
    print("[Paso 4/5] Filtrando y ordenando columnas requeridas...")
    columnas_finales = [
        "califaID",
        "inclination",
        "Re",
        "hi_r",
        "ho_r",
        "r_break",
        "M_bulgue",
        "M_bar",
        "M_disk",
        "pa_calc",
        "h_gas_arcsec",
        "Sigma0_gas",
        "M_gas_1e10",
        "M_star",
        "h_optimo",
    ]

    columnas_faltantes = [col for col in columnas_finales if col not in df_consolidado.columns]
    if columnas_faltantes:
        print(f"[ERROR] Faltan las siguientes columnas requeridas: {columnas_faltantes}")
        sys.exit(1)

    df_final = df_consolidado[columnas_finales].copy()
    df_final = df_final.round(5)

    # Paso 5: Exportación del resultado
    print(f"[Paso 5/5] Exportando tabla consolidada a: {ruta_csv_salida}")
    try:
        ruta_csv_salida.parent.mkdir(parents=True, exist_ok=True)
        df_final.to_csv(ruta_csv_salida, index=False)
        print("  -> Exportación finalizada con éxito.")
    except Exception as e:
        print(f"[ERROR] No se pudo guardar el archivo de salida '{ruta_csv_salida}': {e}")
        sys.exit(1)

    print("\n--- Resumen del DataFrame Final ---")
    print(df_final.to_string(index=False))


if __name__ == "__main__":
    main()
