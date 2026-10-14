#!/usr/bin/env python3
"""
Pipeline Maestro: Reconstrucción completa de los catálogos de datos desde cero.

Orden de ejecución:
  1. src/data_pipeline/filtro_incl.py        -> data/processed/espirales_final.csv
  2. src/data_pipeline/radios.py             -> data/processed/espirales_final_con_bulbo.csv
  3. src/data_pipeline/parametros.py         -> data/processed/parametros.csv
  4. src/data_pipeline/calcular_gas.py       -> data/processed/parametros_con_gas.csv & outputs/figures/perfiles_gas_final.png
  5. src/data_pipeline/calcular_h_optimo.py  -> Cálculo de escalas efectivas h_eff
"""

from pathlib import Path
import sys
import subprocess
import time

BASE_DIR = Path(__file__).resolve().parent

PIPELINE_STEPS = [
    {
        "name": "Paso 1: Selección cinemática y morfológica (filtro_incl.py)",
        "script": BASE_DIR / "src" / "data_pipeline" / "filtro_incl.py",
        "expected_output": BASE_DIR / "data" / "processed" / "espirales_final.csv",
    },
    {
        "name": "Paso 2: Extracción de radios del bulbo (radios.py)",
        "script": BASE_DIR / "src" / "data_pipeline" / "radios.py",
        "expected_output": BASE_DIR / "data" / "processed" / "espirales_final_con_bulbo.csv",
    },
    {
        "name": "Paso 3: Unificación de parámetros y masas estelares (parametros.py)",
        "script": BASE_DIR / "src" / "data_pipeline" / "parametros.py",
        "expected_output": BASE_DIR / "data" / "processed" / "parametros.csv",
    },
    {
        "name": "Paso 4: Inversión de Kennicutt-Schmidt y cálculo de gas (calcular_gas.py)",
        "script": BASE_DIR / "src" / "data_pipeline" / "calcular_gas.py",
        "expected_output": BASE_DIR / "data" / "processed" / "parametros_con_gas.csv",
    },
    {
        "name": "Paso 5: Modelado de escala efectiva de disco (calcular_h_optimo.py)",
        "script": BASE_DIR / "src" / "data_pipeline" / "calcular_h_optimo.py",
        "expected_output": None,
    },
]


def run_pipeline():
    start_total = time.time()
    print("=" * 70)
    print("INICIANDO EJECUCIÓN DEL PIPELINE DE DATOS")
    print(f"Directorio raíz: {BASE_DIR}")
    print("=" * 70)

    for i, step in enumerate(PIPELINE_STEPS, start=1):
        step_name = step["name"]
        script_path = step["script"]
        expected_output = step["expected_output"]

        print(f"\n[{i}/{len(PIPELINE_STEPS)}] {step_name}")
        print(f"  Ejecutando: {script_path.relative_to(BASE_DIR)}")

        if not script_path.exists():
            print(f"  [ERROR] No se encontró el script: {script_path}")
            sys.exit(1)

        t0 = time.time()
        res = subprocess.run([sys.executable, str(script_path)], cwd=str(BASE_DIR))
        dt = time.time() - t0

        if res.returncode != 0:
            print(f"  [FALLO] El script falló con código de salida {res.returncode}")
            sys.exit(res.returncode)

        if expected_output:
            if expected_output.exists():
                print(f"  [OK] Salida verificada: {expected_output.relative_to(BASE_DIR)} ({expected_output.stat().st_size} bytes)")
            else:
                print(f"  [ALERTA] No se encontró el archivo esperado: {expected_output}")
                sys.exit(1)

        print(f"  Completado en {dt:.2f} s")

    total_time = time.time() - start_total
    print("\n" + "=" * 70)
    print(f"PIPELINE COMPLETADO CON ÉXITO en {total_time:.2f} segundos.")
    print("Archivo maestro final generado en: data/processed/parametros_con_gas.csv")
    print("=" * 70)


if __name__ == "__main__":
    run_pipeline()
