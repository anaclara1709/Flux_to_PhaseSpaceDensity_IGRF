import numpy as np
from def_b_igrf import b_igrf
def integrand_k_igrf(s, B_mirror, sol_ode, date, Re):
    """
    A função que será integrada: [B_mirror - B(s)]^(1/2)
    """
    # Obter a posição (X, Y, Z) no ponto 's' através da interpolação densa
    r_at_s = sol_ode.sol(s)
    
    #Conversão ECEF -> CEG para a entrada de b_igrf
    r_mag = np.linalg.norm(r_at_s)
    x_ecef, y_ecef, z_ecef = r_at_s
    lat_rad = np.arcsin(z_ecef / r_mag)
    lon_rad = np.arctan2(y_ecef, x_ecef)

    r, lat_deg, lng_deg = r_mag, np.rad2deg(lat_rad), np.rad2deg(lon_rad)

    #Chamar b_igrf para obter as componentes ECEF
    bx, by, bz = b_igrf(r, lat_deg, lng_deg, date)
    
    # 3. Calcular a Magnitude B(s)
    B_s = np.linalg.norm(np.array([bx, by, bz]))
    
    # Condição para evitar a raiz quadrada de números negativos (B(s) > Bm)
    # Se o ponto de espelho não foi encontrado, isso deve retornar 0.0
    if B_mirror < B_s:
        return 0.0
    else:
        return np.sqrt(B_mirror - B_s)