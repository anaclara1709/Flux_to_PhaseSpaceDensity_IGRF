import numpy as np
from def_b_igrf import b_igrf 
import datetime 
import math

def field_line_igrf_deriv(s, r_vector, B_mirror, date, Re): 
    """
    Função do sistema (derivada) para o algoritmo de integração Runge-Kutta: dr/ds = B/|B|.
    """
    #Extrair as coordenadas de posição do vetor de estado (r_vector). Essa coordenada deve ser cartesiana porque o algoritmo de Runge-Kutta calcula a nova posição usando a derivada cartesiana dr/ds.
    x_ecef, y_ecef, z_ecef = r_vector[0], r_vector[1], r_vector[2]

    #Conversão das coordenadas de ECEF para CEG 
    r_mag = np.linalg.norm(r_vector) 
    lat_rad = np.arcsin(z_ecef / r_mag)
    lon_rad = np.arctan2(y_ecef, x_ecef)

    r = r_mag
    lat = np.rad2deg(lat_rad)
    lng = np.rad2deg(lon_rad)
    
    #Calcular o Campo Magnético  
    bx_ecef, by_ecef, bz_ecef = b_igrf(r, lat, lng, date) 

    #Calcular a magnitude do campo magnetico
    b_vec = np.array([bx_ecef, by_ecef, bz_ecef])
    b_magnitude = np.linalg.norm(b_vec)

    #Calcular a derivada: B / |B| (Vetor Unitário)
    deriv_vector = b_vec / b_magnitude
    
    return np.squeeze(deriv_vector)