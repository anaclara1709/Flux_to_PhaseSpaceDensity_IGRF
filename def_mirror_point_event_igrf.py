import numpy as np
import os
import math
from def_b_igrf import b_igrf

def mirror_point_event_igrf(s, r_vector, B_mirror, date, Re):
    """
    Função de evento que retorna zero quando a Magnitude B (r) alcança B_mirror.
    """
    #Coordenadas em ECEF
    x_ecef, y_ecef, z_ecef = r_vector[0], r_vector[1], r_vector[2]
    
    #Conversão de ECEF para CEG
    r_mag = np.linalg.norm(r_vector) # Isso é 'r' para b_igrf
    lat_rad = np.arcsin(z_ecef / r_mag)
    lon_rad = np.arctan2(y_ecef, x_ecef)

    r = r_mag
    lat = np.rad2deg(lat_rad)
    lng = np.rad2deg(lon_rad)

    #Calcular o Campo Magnético ATUAL B(s) na nova posição
    bx_ecef, by_ecef, bz_ecef = b_igrf(r, lat, lng, date) 
    # Calculo da magnitude do campo 
    b_vec = np.array([bx_ecef, by_ecef, bz_ecef])
    b_magnitude = np.linalg.norm(b_vec)
    #Retornar a diferença: B_mirror - B_current_mag
    # A integração para quando B_current_mag se aproxima de B_mirror
    # Se B_current_mag > B_mirror, o ponto já passou.
    return B_mirror - b_magnitude