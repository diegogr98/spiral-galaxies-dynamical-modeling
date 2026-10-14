"""
Módulo de transformaciones geométricas y proyecciones para dinámica galáctica.
"""

import numpy as np


def rotar_proyeccion(pos, vel, pa_deg, inc_deg):
    """
    Rota posiciones y velocidades al sistema del observador.

    Convención:
      - Disco inicial en plano XY con Lz > 0.
      - Inclinación 'i' rota alrededor del eje X (línea de nodos).
      - PA rota en el plano del cielo alrededor de Z, orientando el eje mayor.

    Parámetros:
        pos (ndarray): Array de posiciones (3, N) en kpc.
        vel (ndarray): Array de velocidades (3, N) en km/s.
        pa_deg (float): Ángulo de posición en grados.
        inc_deg (float): Inclinación en grados.

    Retorna:
        tuple: (pos_p, vel_p) posiciones y velocidades proyectadas (3, N).
    """
    inc = np.radians(inc_deg)
    phi = np.radians(pa_deg - 90.0)

    # 1. Inclinación alrededor de X
    R_x = np.array([
        [1.0, 0.0, 0.0],
        [0.0, np.cos(inc), -np.sin(inc)],
        [0.0, np.sin(inc), np.cos(inc)],
    ])

    # 2. Posición angular en el cielo (alrededor de Z)
    R_z = np.array([
        [np.cos(phi), -np.sin(phi), 0.0],
        [np.sin(phi), np.cos(phi), 0.0],
        [0.0, 0.0, 1.0],
    ])

    R_total = R_z @ R_x
    pos_p = R_total @ pos
    vel_p = R_total @ vel

    return pos_p, vel_p
