import datetime
import numpy as np
import pandas as pd
import os
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from scipy.integrate import quad, solve_ivp
import spacepy.irbempy as ib
import spacepy.time as spt
import spacepy.coordinates as spc
#Somente para o IGRF:
import ppigrf
#Códigos criados: 
from def_plot_cge import plot_cge
from def_ceg_to_ecef import ceg_to_ecef
from plot_3D_field_line import plot_3D_field_line
#Códigos criados somente para o IGRF
from def_b_igrf import b_igrf 
from def_mirror_point_event_igrf import mirror_point_event_igrf
from def_field_line_igrf_deriv import field_line_igrf_deriv
from def_integrand_k_igrf import integrand_k_igrf
plt.ioff()
# -----------------------------------------------------------------------------------------------------------------------
# INPUT
#Raios terrestres: 
Re = 6378
#Informações para rastreio do campo magnetico em coordenadas: 
altitude_re = 1.0 # Altitude desejada acima da superfície
start_latitude_cge = 30.0 # 30.0, exemplo de Latitude Geocêntrica (Equador = 0)
start_longitude_cge = 45.0 # 45.0, exemplo de Longitude Geocêntrica (Greenwich = 0)
start_altitude_km_cge = Re + (altitude_re * Re)
#Lista de pitch angles
pitch_angles_deg = np.arange(0.5, 90.1, 0.5)
#K fixo desejado
K_fixo_desejado_nThalf_km = 1475038.5596173669
#Ano da análise
field_year = 2015
#Data de início e fim da análise (o período de análise deve ser de um dia)
start_date = datetime.datetime(field_year, 1, 1, 0, 0, 0)
end_date = datetime.datetime(field_year, 1, 2, 0, 0, 0)
#Intervalo de tempo 
time_delta = datetime.timedelta(hours=1)
#Diretorio base para armazenar o resultado da análise
day4path = start_date.strftime('%Y%m%d') #dia para armazenar 
base_dir = f"D:/INPE_PCI/SIMULATION OUTER RADIATION BELT 1D/initial condition/github publico (igrf)/results/{field_year}/{day4path}"
# -----------------------------------------------------------------------------------------------------------------------
#Criar todos os diretorios para armazenar resultados da fase 1
directory_path = os.path.join(base_dir, 'STEP 1')
os.makedirs(directory_path, exist_ok=True)
print(f"Verificando/Criando a pasta: '{directory_path}'")
directory_path_plots = os.path.join(directory_path, 'PLOTS')
os.makedirs(directory_path_plots, exist_ok=True)
print(f"Verificando/Criando a pasta: '{directory_path_plots}'")
#Criar diretorio para salvar testes
directory_tests_path = os.path.join(base_dir, 'STEP 1_TESTS')
os.makedirs(directory_tests_path, exist_ok=True)
# -----------------------------------------------------------------------------------------------------------------------
#Figura para testar a localização inicial em coordenadas cge
image_start_cge = "start_cge.png"
path_start_cge = os.path.join(directory_tests_path, image_start_cge)
print(f"Gerando plot para r = {start_altitude_km_cge/Re:.1f} RE...")
plot_cge(Re=Re, r=start_altitude_km_cge, theta_deg=start_latitude_cge, phi_deg=start_longitude_cge, save_path=path_start_cge)
print(f"Processo de plotagem concluído. Verificar imagem na pasta: {directory_tests_path}")
# -----------------------------------------------------------------------------------------------------------------------
#Conversão das coordenadas de ceg para ecef
x_ecef, y_ecef, z_ecef = ceg_to_ecef(r=start_altitude_km_cge, theta_deg=start_latitude_cge, phi_deg=start_longitude_cge)
#PARÂMETROS PARA O CALCULO DA LINHA DE CAMPO:
#O s_max deve ser grande o suficiente para atingir o ponto de espelho. Um valor seguro é 5-10 RE, mas usaremos 100.000 km como um valor grande.
s_max = 5000000.0  # 500k km é uma distância segura
rtol=1e-5
atol=1e-8
# -----------------------------------------------------------------------------------------------------------------------
# Declarar dicionario com os resultados de K para cada instante do intervalo temporal
k_results = {} 
# -----------------------------------------------------------------------------------------------------------------------
# Gerar uma lista de todos os instantes de tempo a serem analisados
current_date = start_date
dates_to_analyze = []
while current_date <= end_date:
    dates_to_analyze.append(current_date)
    current_date += time_delta
# -----------------------------------------------------------------------------------------------------------------------
#LOOP DATE
for date in dates_to_analyze:
    print(f"\n=======================================================")
    print(f"PROCESSANDO INSTANTE DE TEMPO: {date.strftime('%Y-%m-%d %H UT')}")
    print(f"=======================================================")
    #Variavel date em string para criar a chave do dicionário k_results
    date4loop=date.strftime('%Y%m%d%H%M%S')
    #Calculo o b usando o igrf 
    bx_ecef, by_ecef, bz_ecef = b_igrf(start_altitude_km_cge , start_latitude_cge, start_longitude_cge, date)
    #Calculo da magnitude 
    b_magnitude = math.sqrt((bx_ecef ** 2) + (by_ecef ** 2) + (bz_ecef ** 2))
    #Posição inicial (vetor [X, Y, Z]) em km
    r_start_vector = np.array([x_ecef, y_ecef, z_ecef])
    #Declarar a lista de K para cada pitch angle 
    k_values = []
    #---- LOOP PITCH ANGLES
    for alpha_deg in pitch_angles_deg:
        #Pitch angle em radianos
        alpha_rad = np.deg2rad(alpha_deg)
        # Calcular o campo de espelho Bm 
        B_mirror = b_magnitude / (np.sin(alpha_rad)**2)
        # Configurações adicionais para o evento:
        mirror_point_event_igrf.terminal = True # O integrador para quando o evento ocorre
        mirror_point_event_igrf.direction = -1  # A integração para quando a função vai de >0 para <0 (B_atual excede B_m)
        # --------------
        #RASTREAMENTO -- Argumentos adicionais para a função de derivada
        deriv_args = (B_mirror, date, Re)
        #Rastreamento para o NORTE: 
        s_span_norte = [0, s_max] #Intervalo de integração
        print("\nIniciando Rastreamento da Linha de Campo (Direção Norte) com evento de Parada...")
            #Runge-Kutta NORTE:
        sol_norte = solve_ivp(field_line_igrf_deriv, s_span_norte, r_start_vector, args=deriv_args, method='RK45', dense_output=True, rtol=rtol, atol=atol, events=mirror_point_event_igrf)
        print(f"Rastreamento Finalizado. Ponto de Espelho NORTE encontrado em s = {sol_norte.t_events[0][0]:.2f} km")
        #Rastreamento para o SUL: 
        s_span_sul = [0, -s_max]
        print("\nIniciando Rastreamento da Linha de Campo (Direção Sul) com evento de Parada...")
            #Runge-Kutta SUL
        sol_sul = solve_ivp(field_line_igrf_deriv, s_span_sul, r_start_vector, args=deriv_args, method='RK45', dense_output=True, rtol=rtol, atol=atol, events=mirror_point_event_igrf)
        if sol_sul.t_events[0].size > 0:
            print(f"Rastreamento Finalizado. Ponto de Espelho SUL encontrado em s = {sol_sul.t_events[0][0]:.2f} km")
        else:
            s_m_sul = 0.0 # Se o evento não for encontrado
            print("Aviso: Ponto de Espelho Sul não foi alcançado no limite de integração.")
        # --------------
        #PLOT 3D DE LINHAS DE CAMPO 
        #Perspectivas que serão salvas
        perspectivas = [
            {'elev': 20, 'azim': 45},   # Vista padrão
            {'elev': 90, 'azim': 0},    # Vista do Polo (Topo)
            {'elev': 0, 'azim': 0},     # Vista do Plano XZ (Meridiano de Greenwich)
            {'elev': 0, 'azim': 90}     # Vista do Plano YZ (Plano Equatorial - 90 E)
        ]
        print("\nIniciando Geração de plots 3D das linhas de campo de múltiplas Perspectivas...")
        for p in perspectivas:
            plot_3D_field_line(sol_norte=sol_norte, sol_sul=sol_sul, r_start_vector=r_start_vector, Re=Re, alpha_deg=alpha_deg, elev=p['elev'], azim=p['azim'], save_dir=directory_tests_path, date=date)
        print("\nTodas as perspectivas foram salvas no diretório de testes.")
        #----------------
        #CALCULO DE K
        #Extração dos Limites de Integração
        s_m_norte = sol_norte.t_events[0][0] if sol_norte.t_events[0].size > 0 else 0.0
        s_m_sul = sol_sul.t_events[0][0] if sol_sul.t_events[0].size > 0 else 0.0
        if s_m_norte != 0.0 and s_m_sul != 0.0:
            # A linha está fechada e o K pode ser calculado        
            # I. Integral Norte
            integral_norte, erro_norte = quad(integrand_k_igrf, 0, s_m_norte, args=(B_mirror, sol_norte, date, Re))
            # II. Integral Sul
            integral_sul, erro_sul = quad(integrand_k_igrf, s_m_sul, 0, args=(B_mirror, sol_sul, date, Re))
            #Calculo de K_total 
            K_total = integral_norte + integral_sul
        else:
            # Caso a linha esteja aberta ou tenha falhado
            K_total = np.nan
            print(f"AVISO: Linha aberta ou falha no rastreamento para alpha={alpha_deg}")
        k_values.append(K_total)
    # ---- FIM DO LOOP PITCH ANGLES
    #----------------
    #Converter para DataFrame
    df_k_alpha = pd.DataFrame({'Alpha_deg': pitch_angles_deg, 'K_invariante': k_values})
    #Filtrar valores válidos (onde o K foi calculado, K é real)
    df_valid = df_k_alpha.dropna(subset=['K_invariante'])
    if df_valid.shape[0] < 2:
        print(f"\nERRO: Poucos pontos (N={df_valid.shape[0]}) para interpolação K(alpha).")
    else:
        #Criar a função de interpolação: alpha = f(K)
        #Interpolação da dependência: Alpha (eixo Y) em função de K (eixo X)
        f_alpha_of_k = interp1d(df_valid['K_invariante'], df_valid['Alpha_deg'], kind='linear', fill_value='extrapolate')
        # Determinar o ângulo de pitch fixo (alpha_K)
        K_target = float(K_fixo_desejado_nThalf_km)
        alpha_k_result = f_alpha_of_k(K_target)
        print(f"\n--- Resultado final (interpolação) ---")
        print(f"K fixo alvo: {K_target:.2e} (nT^1/2 * km)")
        #----------------
        #CRIAR PLOT DE VALIDAÇÃO
        plt.figure(figsize=(8, 6))    
        # Plotar os Pontos Calculados K(alpha)
        plt.plot(df_valid['Alpha_deg'], df_valid['K_invariante'], 'k-o', label=r'K($\alpha$) Calculado', linewidth=1, markersize=3)
        # Plotar o Ponto de Interpolação (Alpha_K Fixo)
        plt.plot(alpha_k_result, K_target, 'ro', label=r'$\alpha_K$ Fixo ({:.2f} deg)'.format(alpha_k_result), markersize=8)
        # Adicionar Linhas de Referência
            # Linha vertical 
        plt.axvline(x=alpha_k_result, color='r', linestyle='--', linewidth=0.8, alpha=0.6)
            # Linha horizontal 
        plt.axhline(y=K_target, color='r', linestyle='--', linewidth=0.8, alpha=0.6)
        #Configuração do Gráfico
        plt.xlabel(r'Pitch angle ($\alpha$) (degree)')
        plt.ylabel(r'Second Adiabatic Invariant $K$ ($nT^{1/2} \cdot km$)')
        plt.title(f'Function K($\\alpha$) and Interpolation (Time: {date.strftime("%Y-%m-%d %H UT")})')
        plt.grid(True, linestyle=':', alpha=0.7)
        plt.legend()
        #Salvar e Fechar
        plot_k_alpha_filename = f"K_alpha_curve_{date4loop}_Kfixo_{str(K_fixo_desejado_nThalf_km)}.png"
        plot_k_alpha_path = os.path.join(directory_path_plots, plot_k_alpha_filename)
        plt.savefig(plot_k_alpha_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Curva K($\\alpha$) salva em: {plot_k_alpha_path}")
        #----------------
        #CRIAR PLOT DE VALIDAÇÃO (SQRT(G) * R_E)
        # Fator de conversão: (sqrt(nT) / sqrt(G)) * (1/km)
        # sqrt(G) = sqrt(10^5 * nT) => sqrt(nT)/sqrt(G) = 1 / sqrt(10^5) approx 0.00316
        conversion_factor = 1.0 / np.sqrt(1e5) # Fator para converter nT^1/2 -> G^1/2
        # Conversão das colunas válidas
        df_valid['K_RE_G'] = df_valid['K_invariante'] / Re * conversion_factor
        K_target_RE_G = K_target / Re * conversion_factor
        alpha_k_result_re_g = f_alpha_of_k(K_target) # O alpha_k_result é o mesmo
        plt.figure(figsize=(8, 6))
        # Plotar os Pontos Calculados K(alpha)
        # X = Alpha, Y = K_RE_G (Novo Eixo)
        plt.plot(df_valid['Alpha_deg'], df_valid['K_RE_G'], 'k-o', label=r'K($\alpha$) Calculado', linewidth=1, markersize=3)
        # Plotar o Ponto de Interpolação (Alpha_K Fixo)
        # X = alpha_k_result, Y = K_target_RE_G
        plt.plot(alpha_k_result_re_g, K_target_RE_G, 'ro', label=r'$\alpha_K$ Fixo ({:.2f} deg)'.format(alpha_k_result_re_g), markersize=8)
        # Adicionar Linhas de Referência
        plt.axvline(x=alpha_k_result_re_g, color='r', linestyle='--', linewidth=0.8, alpha=0.6)
        plt.axhline(y=K_target_RE_G, color='r', linestyle='--', linewidth=0.8, alpha=0.6)
        # Configuração do Gráfico
        plt.xlabel(r'Pitch angle ($\alpha$) (degree)')
        # Novo rotulo de unidade
        plt.ylabel(r'Second Adiabatic Invariant $K$ ($\sqrt{G} \cdot R_E$)') 
        plt.title(f'Função K($\\alpha$) e Interpolação (Unidade Clássica - {date.strftime("%Y-%m-%d %H UT")})')
        plt.grid(True, linestyle=':', alpha=0.7)
        plt.legend()
        # Salvar e Fechar
        plot_k_alpha_filename_G = f"K_alpha_curve_{date4loop}_Kfixo_{str(K_fixo_desejado_nThalf_km)}_sqrt(G)Re.png"
        plot_k_alpha_path_G = os.path.join(directory_path_plots, plot_k_alpha_filename_G)
        plt.savefig(plot_k_alpha_path_G, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Curva K($\\alpha$) (Unidade G) salva em: {plot_k_alpha_path_G}")
        #----------------
    #Armazenar o resultado final do dia no dicionário (Chave = Data)
    k_results[date4loop] = {'alpha_K': alpha_k_result, 'K_input': K_fixo_desejado_nThalf_km}

#----------------
#SALVAR CSV
print("\nSalvar resultados finais da série temporal...")
df_k_results_final = pd.DataFrame.from_dict(k_results, orient='index')
df_k_results_final.index.name = 'Data'

csv_filename = f"STEP_1_{day4path}_Kfixo_{str(K_fixo_desejado_nThalf_km)}.csv"
csv_path = os.path.join(directory_path, csv_filename)
df_k_results_final.to_csv(csv_path)
print(f"Resultados K salvos com sucesso em: {csv_path}")












exit()
#Calculo o b usando o igrf (VARIAR A DATA)
bx_ecef, by_ecef, bz_ecef = b_igrf(start_altitude_km_cge , start_latitude_cge, start_longitude_cge, end_date)
print(bx_ecef, by_ecef, bz_ecef)
#Calculo da magnitude 
b_magnitude = math.sqrt((bx_ecef ** 2) + (by_ecef ** 2) + (bz_ecef ** 2))
print(b_magnitude)
# -----------------------------------------------------------------------------------------------------------------------

#1) Pontos espelho:
# O ângulo de pitch alfa deve ser testado. Vamos escolher um para o teste:
alpha_deg_teste = 45.0
alpha_rad_teste = np.deg2rad(alpha_deg_teste)
# Calcular o campo de espelho Bm 
B_mirror = b_magnitude / (np.sin(alpha_rad_teste)**2)
print(f"Campo de Espelho (Bm) para alpha={alpha_deg_teste}°: {B_mirror:.2f} nT")
# Configurações adicionais para o evento:
mirror_point_event_igrf.terminal = True # O integrador para quando o evento ocorre
mirror_point_event_igrf.direction = -1  # A integração para quando a função vai de >0 para <0 (B_atual excede B_m)

#2)Rastreamento para o Norte: 
    #a) Posição inicial (vetor [X, Y, Z]) em km
r_start_vector = np.array([x_ecef, y_ecef, z_ecef])
    #b) Intervalo de integração: [0, s_max]. O s_max deve ser grande o suficiente para atingir o ponto de espelho. Um valor seguro é 5-10 RE, mas usaremos 100.000 km como um valor grande.
s_max = 5000000.0 
s_span_norte = [0, s_max]
    #c) Argumentos adicionais para a função de derivada
deriv_args = (B_mirror, end_date, Re)
print("\nIniciando Rastreamento da Linha de Campo (Direção Norte) com evento de Parada...")
    #d) Runge-Kutta 
sol_norte = solve_ivp(
    field_line_igrf_deriv, 
    s_span_norte, 
    r_start_vector, 
    args=deriv_args, 
    method='RK45', 
    dense_output=True,
    rtol=1e-5, 
    atol=1e-8, 
    events=mirror_point_event_igrf
)
# A linha de campo rastreada está em sol_norte_final.y
print(f"Rastreamento Finalizado. Ponto de Espelho NORTE encontrado em s = {sol_norte.t_events[0][0]:.2f} km")
#3) Rastreamento para o Sul:
s_span_sul = [0, -s_max]
print("\nIniciando Rastreamento da Linha de Campo (Direção Sul) com evento de Parada...")
sol_sul = solve_ivp(
    field_line_igrf_deriv, 
    s_span_sul, 
    r_start_vector, 
    args=deriv_args, 
    method='RK45', 
    dense_output=True,
    rtol=1e-5, 
    atol=1e-8,
    events=mirror_point_event_igrf
)
if sol_sul.t_events[0].size > 0:
    print(f"Rastreamento Finalizado. Ponto de Espelho SUL encontrado em s = {sol_sul.t_events[0][0]:.2f} km")
else:
    s_m_sul = 0.0 # Se o evento não for encontrado
    print("Aviso: Ponto de Espelho Sul não foi alcançado no limite de integração.")

#----------------
#PLOT 3D DE LINHAS DE CAMPO MAGNÉTICO
#Perspectivas que serão salvas
perspectivas = [
    {'elev': 20, 'azim': 45},   # Vista padrão
    {'elev': 90, 'azim': 0},    # Vista do Polo (Topo)
    {'elev': 0, 'azim': 0},     # Vista do Plano XZ (Meridiano de Greenwich)
    {'elev': 0, 'azim': 90}     # Vista do Plano YZ (Plano Equatorial - 90 E)
]
print("\nIniciando Geração de plots 3D das linhas de campo de múltiplas Perspectivas...")
for p in perspectivas:
    plot_3D_field_line(
        sol_norte=sol_norte, 
        sol_sul=sol_sul, 
        r_start_vector=r_start_vector, 
        Re=Re, 
        alpha_deg=alpha_deg_teste, 
        elev=p['elev'], 
        azim=p['azim'],
        save_dir=directory_tests_path # Pasta de testes
    )
print("\nTodas as perspectivas foram salvas no diretório de testes.")
#----------------
#Calculo de K
#Extração dos Limites de Integração
s_m_norte = sol_norte.t_events[0][0]
s_m_sul = sol_sul.t_events[0][0] 
#Executar Quadratura (Integral de s'_m a s_m)
# A integração é feita em duas partes se os pontos forem simétricos, mas o scipy.integrate.quad pode integrar diretamente no intervalo [s_m_sul, s_m_norte].
# A integral total é dividida em duas chamadas para garantir que o vetor sol_ode seja correto (sol_norte é válido para s>=0, sol_sul é válido para s<=0).
# I. Integral Norte (0 até s_m)
integral_norte, erro_norte = quad(integrand_k_igrf, 0, s_m_norte, args=(B_mirror, sol_norte, end_date, Re))
# II. Integral Sul (s'_m até 0)
# A integral é feita no valor absoluto da distância, mas o intervalo deve ser [s_m_sul, 0]
integral_sul, erro_sul = quad(integrand_k_igrf, s_m_sul, 0, args=(B_mirror, sol_sul, end_date, Re))
#Calculo de K_total 
K_total = integral_norte + integral_sul
print(K_total)