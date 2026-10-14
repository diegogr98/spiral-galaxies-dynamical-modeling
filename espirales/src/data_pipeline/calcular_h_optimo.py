"""
Módulo del pipeline: Cálculo de escala de longitud efectiva (h_eff) para galaxias con quiebre de disco (Type II/III).
Basado en el formalismo de Méndez-Abreu et al.
"""

from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
PROCESSED_DIR = BASE_DIR / "data" / "processed"


def calculate_effective_h(h_in, h_out, r_break, r_max=40.0, num_points=1000):
    """
    Evalúa la ecuación de Méndez-Abreu et al. y ajusta un único disco exponencial.
    """
    r = np.linspace(0.1, r_max, num_points)

    # Factor de continuidad en r_break
    continuity_factor = np.exp(-r_break * (h_out - h_in) / (h_out * h_in))

    # Ecuación del perfil
    theta = (r < r_break).astype(float)
    I_profile = np.exp(-r / h_in) * theta + continuity_factor * np.exp(-r / h_out) * (1.0 - theta)

    # Ponderación por área del disco (pesos proporcionales a r)
    weights = np.sqrt(r)

    # Ajuste lineal en espacio logarítmico: ln(I) = b + m*r, donde m = -1 / h_eff
    fit_params = np.polyfit(r, np.log(I_profile), deg=1, w=weights)
    m = fit_params[0]
    h_eff = -1.0 / m

    # Ajuste sin pesos para comparar
    fit_simple = np.polyfit(r, np.log(I_profile), deg=1)
    h_eff_simple = -1.0 / fit_simple[0]

    return h_eff, h_eff_simple


def main():
    print("--- Ejecutando calcular_h_optimo.py ---")
    galaxies = {
        684: {"hi": 19.54466, "ho": 9.42834, "r_break": 30.50633, "rmax": 45.0},
        813: {"hi": 26.50555, "ho": 25.08450, "r_break": 113.0, "rmax": 50.0},
    }

    # Si existe el archivo parametros.csv, podemos verificar valores
    path_param = PROCESSED_DIR / "parametros.csv"
    if path_param.exists():
        df = pd.read_csv(path_param)
        for cid in list(galaxies.keys()):
            sub = df[df["califaID"] == cid]
            if not sub.empty:
                row = sub.iloc[0]
                if pd.notna(row.get("hi_r")) and pd.notna(row.get("ho_r")) and pd.notna(row.get("r_break")):
                    galaxies[cid]["hi"] = float(row["hi_r"])
                    galaxies[cid]["ho"] = float(row["ho_r"])
                    galaxies[cid]["r_break"] = float(row["r_break"])

    print("--- RESULTADOS DE ESCALA EFECTIVA ÚNICA (h_eff) ---")
    for cid, params in galaxies.items():
        h_eff, h_simple = calculate_effective_h(
            params["hi"], params["ho"], params["r_break"], r_max=params["rmax"]
        )
        print(f"CALIFA {cid}:")
        print(
            f"  h_in = {params['hi']:.2f}'', h_out = {params['ho']:.2f}'', r_break = {params['r_break']:.1f}''"
        )
        print(f"  -> h_eff (ponderado por área) = {h_eff:.2f}''")
        print(f"  -> h_eff (ajuste lineal simple) = {h_simple:.2f}''\n")


if __name__ == "__main__":
    main()