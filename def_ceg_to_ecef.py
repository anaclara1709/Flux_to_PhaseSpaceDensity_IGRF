import numpy as np

def ceg_to_ecef(r, theta_deg, phi_deg):
    """
    Converte coordenadas Esféricas Geocêntricas (CEG) para ECEF (Cartesianas).
    
    CEG:
    - r: Raio Geocêntrico (distância radial em km). Esse valor deve ser Re + (altitude * Re).  
    - theta_deg: Latitude Geocêntrica (em graus, ângulo a partir do Equador).
    - phi_deg: Longitude Geocêntrica (em graus, ângulo a partir de Greenwich).
    
    ECEF: (X, Y, Z) Cartesianas.

    Args:
        r (float): Raio Geocêntrico (unidade de distância em km).
        theta_deg (float): Latitude Geocêntrica (em graus).
        phi_deg (float): Longitude Geocêntrica (em graus).

    Returns:
        tuple: Coordenadas ECEF Cartesianas (X, Y, Z).
    """
    # 1. Converter ângulos de graus para radianos
    theta_rad = np.deg2rad(theta_deg)
    phi_rad = np.deg2rad(phi_deg)
    
    # 2. Aplicar as fórmulas de conversão
    # X = r * cos(theta) * cos(phi)  (Alinhado com Greenwich e Equador)
    X = r * np.cos(theta_rad) * np.cos(phi_rad)
    
    # Y = r * cos(theta) * sin(phi)  (Alinhado com 90° Leste)
    Y = r * np.cos(theta_rad) * np.sin(phi_rad)
    
    # Z = r * sin(theta)          (Alinhado com o Polo Norte)
    Z = r * np.sin(theta_rad)
    
    return X, Y, Z