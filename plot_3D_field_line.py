import datetime
import os
import numpy as np
import matplotlib.pyplot as plt

def plot_3D_field_line(sol_norte, sol_sul, r_start_vector, Re, alpha_deg, elev, azim, save_dir,date):
    """
    Gera e salva o gráfico da linha de campo em uma perspectiva específica.
    """

    
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # -----------------------------------------------------------------------------------------
    # PREPARAÇÃO DE DADOS 
    # -----------------------------------------------------------------------------------------
    # Certifique-se de que todas as variáveis necessárias (Re, r_start_vector, etc.)
    # são passadas para esta função.
    
    x_sul = sol_sul.y[0][::-1]
    y_sul = sol_sul.y[1][::-1]
    z_sul = sol_sul.y[2][::-1]
    
    x_norte = sol_norte.y[0]
    y_norte = sol_norte.y[1]
    z_norte = sol_norte.y[2]
    
    X_start, Y_start, Z_start = r_start_vector
    r_start_re = np.linalg.norm(r_start_vector) / Re
    max_len = max(np.linalg.norm(r_start_vector), 2 * Re) * 1.5
    
    # -----------------------------------------------------------------------------------------
    # PLOTAGEM (Esfera, Eixos, Linhas, etc.) 
    # -----------------------------------------------------------------------------------------
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_terra = Re * np.outer(np.cos(u), np.sin(v))
    y_terra = Re * np.outer(np.sin(u), np.sin(v))
    z_terra = Re * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x_terra, y_terra, z_terra, color='skyblue', alpha=0.3, linewidth=0)

    # Eixos Centrais 
    ax.plot([-max_len, max_len], [0, 0], [0, 0], color='gray', linewidth=0.8, linestyle='-', alpha=0.7)
    ax.plot([0, 0], [-max_len, max_len], [0, 0], color='gray', linewidth=0.8, linestyle='-', alpha=0.7)
    ax.plot([0, 0], [0, 0], [-max_len, max_len], color='gray', linewidth=0.8, linestyle='-', alpha=0.7)

    # Linhas de Campo
    ax.plot(x_norte, y_norte, z_norte, color='blue', linewidth=2, label='Linha de Campo (Norte)')
    ax.plot(x_sul, y_sul, z_sul, color='green', linewidth=2, label='Linha de Campo (Sul)')
    ax.scatter(X_start, Y_start, Z_start, color='red', s=50, label=f'Ponto Inicial (r={r_start_re:.1f} $R_E$)', marker='o')

    # ANOTAÇÕES E CONFIGURAÇÃO
    label_offset = max_len * 0.98 
    ax.text(label_offset, 0, 0, 'X ($R_E$)', color='black', fontsize=12, horizontalalignment='left')
    ax.text(0, label_offset, 0, 'Y ($R_E$)', color='black', fontsize=12, horizontalalignment='left')
    ax.text(0, 0, label_offset, 'Z ($R_E$)', color='black', fontsize=12, horizontalalignment='center')
        
    ax.set_axis_off() 
    ax.set_xlim([-max_len, max_len])
    ax.set_ylim([-max_len, max_len])
    ax.set_zlim([-max_len, max_len])
    ax.set_box_aspect((1, 1, 1))
    ax.set_title(f'Linha de Campo Magnético (TS05) - $\\alpha_0$={alpha_deg}°', fontsize=14)
    ax.legend(loc='upper right')
    
    # VARIAÇÃO DA PERSPECTIVA E SALVAMENTO MÚLTIPLO
    ax.view_init(elev=elev, azim=azim) 
    
    # Cria o nome do arquivo incluindo os ângulos
    day4path = date.strftime('%Y%m%d%H%M%S')
    filename = f"Field_line_{day4path}_PitchAngle_{int(alpha_deg)}_elev{elev}_azim{azim}.png"
    full_path = os.path.join(save_dir, filename)
    
    # SALVA
    plt.savefig(full_path, dpi=300, bbox_inches='tight')
    print(f"Salvo: {filename}")
    
    # Fecha a figura para liberar memória
    plt.close(fig)