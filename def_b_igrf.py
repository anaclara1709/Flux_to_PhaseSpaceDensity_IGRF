import numpy as np
import ppigrf
import spacepy.irbempy as ib
import spacepy.time as spt
import spacepy.coordinates as spc
def b_igrf(r, lat, lng, date): 
    """
    Calcula o campo magnetico usando o igrf.
    Args:
        r (float), lat (deg) e lng: coordenada CEG (Coordenada Esferica Geocentrica).
        date: A entrada deve estar no formato datetime.datetime(ANO, MES, DIA, HORA, MINUTO, SEGUNDO). 

    Return:
        Campo magnetico em coordendas ECEF. 
    """
    #colatitude
    co_lat = 90.0 - lat
    #calculo do campo magnetico usando o igrf (coordenadas geocentricas)
    Br, Btheta, Bphi = ppigrf.igrf_gc(r=np.array([r]), theta=np.array([co_lat]), phi=np.array([lng]), date=date)
    Br_val = Br[0]
    Btheta_val = Btheta[0]
    Bphi_val = Bphi[0]
    #converter vetor CEG (Br, Btheta, Bphi) para ECEF (Bx, By, Bz)  por meio da matriz de rotação é Br, Btheta (polar), Bphi.
    #Projeção Radial + Projeção Polar + Projeção Azimutal
    # A projeção B_r * sin(theta) é a componente radial na direção Z.
    # A projeção B_theta * cos(theta) é a componente polar na direção Z.
    # Como Btheta é Sul, ele contribui NEGATIVAMENTE para o Eixo Z (Norte).
    theta_rad = np.deg2rad(lat) 
    phi_rad = np.deg2rad(lng)
    # Bx_ecef
    bx_ecef = Br_val * np.cos(theta_rad) * np.cos(phi_rad) + Btheta_val * np.sin(theta_rad) * np.cos(phi_rad) - Bphi_val * np.sin(phi_rad)
    # By_ecef
    by_ecef = Br_val * np.cos(theta_rad) * np.sin(phi_rad) + Btheta_val * np.sin(theta_rad) * np.sin(phi_rad) + Bphi_val * np.cos(phi_rad)
    # Bz_ecef
    bz_ecef = Br_val * np.sin(theta_rad) - Btheta_val * np.cos(theta_rad)

    return bx_ecef, by_ecef, bz_ecef