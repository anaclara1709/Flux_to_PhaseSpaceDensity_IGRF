import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_cge(Re, r, theta_deg, phi_deg, save_path):
    """
    #Código que gera figura para representar um ponto em coordenadas esféricas geocêntricas 

    #Args:
    #Re (float): raio terrestre (em km). 
    #r (float): Raio Geocêntrico (unidade de distância em Re). Esse valor deve ser Re + (altitude * Re).  
    #theta_deg (float): Latitude Geocêntrica (em graus, 0 no Equador). Unidade de entrada graus.
    #phi_deg (float): Longitude Geocêntrica (em graus, 0 em Greenwich). Unidade de entrada graus.
    """
# Converter graus para radianos
    theta_rad = np.deg2rad(theta_deg)
    phi_rad = np.deg2rad(phi_deg)
    # Fórmulas de conversão (usando theta = Latitude)
    X = r * np.cos(theta_rad) * np.cos(phi_rad)
    Y = r * np.cos(theta_rad) * np.sin(phi_rad)
    Z = r * np.sin(theta_rad)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
        
    # 1. DEFINIÇÕES DE ESCALA E COMPRIMENTO
    max_val_plot = max(r, Re * 2) # Garante que pelo menos 2 RE sejam visíveis
    max_len = max_val_plot * 1.5 # Aumenta o comprimento total da linha do eixo
    # Gerar pontos de tick em múltiplos de RE (1 RE, 2 RE, 3 RE, etc.)
    num_re_max = int(np.ceil(max_val_plot / Re))
    re_ticks_vals = np.arange(1, num_re_max + 1) * Re# [1*RE, 2*RE, 3*RE, ...]
        
    # Adicionar o tick do ponto, se ele não for um RE exato ou se for maior
    if r > Re * num_re_max:
            re_ticks_vals = np.sort(np.unique(np.append(re_ticks_vals, r)))
        
    # Offset para posicionar os rótulos dos números dos REs
    text_offset_re = max_val_plot * 0.05  # Ajustado com base em max_val_plot
        
    # 2. DESENHAR OS EIXOS CENTRAIS MANUAIS (cinza para o fundo)
    #X Axis
    ax.plot([-max_len, max_len], [0, 0], [0, 0], color='gray', linewidth=0.8, linestyle='-', alpha=0.7)
    # Y Axis
    ax.plot([0, 0], [-max_len, max_len], [0, 0], color='gray', linewidth=0.8, linestyle='-', alpha=0.7)
    # Z Axis
    ax.plot([0, 0], [0, 0], [-max_len, max_len], color='gray', linewidth=0.8, linestyle='-', alpha=0.7)
        
    # 3. Desenhar a Esfera da Terra (Raio = Re)
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_terra = Re * np.outer(np.cos(u), np.sin(v))
    y_terra = Re * np.outer(np.sin(u), np.sin(v))
    z_terra = Re * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x_terra, y_terra, z_terra, color='skyblue', alpha=0.4, linewidth=0)
        
    # 4. PLOTAR O VETOR RAIO
    ax.plot([0, X], [0, Y], [0, Z], color='black', linestyle='-', linewidth=1.5, alpha=0.9, label='Vetor Raio (r)')

    # 5. Plotar o Ponto
    ax.scatter(X, Y, Z, color='black', s=30, label=f'Ponto (r={r/Re:.1f} $R_E$)', marker='o')

    # 6. ANOTAÇÕES NOS EIXOS 
    # Coordenada um pouco além do final do eixo para o rótulo
    label_offset = max_len * 0.95    
    # Rótulos nas extremidades 
    # Eixo X
    ax.text(label_offset, 0, 0, 'X ($R_E$)', color='black', fontsize=12, horizontalalignment='left')
    # Eixo Y
    ax.text(0, label_offset, 0, 'Y ($R_E$)', color='black', fontsize=12, horizontalalignment='left')
    # Eixo Z
    ax.text(0, 0, label_offset, 'Z ($R_E$)', color='black', fontsize=12, horizontalalignment='center')
        
    # Marcar os raios terrestres (Solicitação 1)
    for tick_val in re_ticks_vals:
    # Apenas se o tick for maior que 0
        if tick_val > 0:
            re_label = f'{tick_val/Re:.0f}' if np.isclose(tick_val/Re, np.round(tick_val/Re)) else f'{tick_val/Re:.1f}'
                
            # Marcador no Eixo X
            ax.text(tick_val, 0, 0, re_label, color='black', fontsize=10, horizontalalignment='center', verticalalignment='top')
            ax.text(-tick_val, 0, 0, f'-{re_label}', color='black', fontsize=10, horizontalalignment='center', verticalalignment='top')
                
            # Marcador no Eixo Y
            ax.text(0, tick_val, 0, re_label, color='black', fontsize=10, horizontalalignment='center', verticalalignment='top')
            ax.text(0, -tick_val, 0, f'-{re_label}', color='black', fontsize=10, horizontalalignment='center', verticalalignment='top')
                
            # Marcador no Eixo Z
            ax.text(0, 0, tick_val, re_label, color='black', fontsize=10, horizontalalignment='center', verticalalignment='bottom')
            ax.text(0, 0, -tick_val, f'-{re_label}', color='black', fontsize=10, horizontalalignment='center', verticalalignment='top')

    # 7. CONFIGURAÇÃO FINAL

    # Desabilitar a caixa, grade e ticks padrão
    ax.set_axis_off() 
        
    # Definir limites e garantir o aspecto igual
    ax.set_xlim([-max_len, max_len])
    ax.set_ylim([-max_len, max_len])
    ax.set_zlim([-max_len, max_len])
        
    ax.set_box_aspect((1, 1, 1))
        
    ax.set_title('Visualização do Ponto no Sistema Geocêntrico (Esquemático)', fontsize=14)
    ax.legend(loc='upper right')
        
    # Ajustar o ângulo de visão para se aproximar do seu esboço
    ax.view_init(elev=30, azim=-60) 
    
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Plot salvo em: {save_path}")
        
    plt.close(fig)
