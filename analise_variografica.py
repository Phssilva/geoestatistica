import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import curve_fit
import os

# Configurações para os gráficos
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Criar diretório para salvar as figuras se não existir
if not os.path.exists('figuras'):
    os.makedirs('figuras')

# Carregar os dados limpos
df = pd.read_csv('dados_potassio_limpos.csv')

# Extrair coordenadas e valores
coords = df[['x', 'y']].values
values = df['K_mg'].values

# Função para calcular o variograma experimental
def calcular_variograma_experimental(coords, values, n_lags=15, max_dist=None, tolerance=0.5, angle=0, angle_tolerance=45):
    """
    Calcula o variograma experimental.
    
    Parâmetros:
    - coords: array com coordenadas [x, y]
    - values: array com valores da variável
    - n_lags: número de passos (lags)
    - max_dist: distância máxima a considerar
    - tolerance: tolerância do passo (lag tolerance)
    - angle: ângulo da direção (em graus, 0 = leste, 90 = norte)
    - angle_tolerance: tolerância angular (em graus)
    
    Retorna:
    - lag_distances: distâncias médias para cada lag
    - gamma: valores do semivariograma para cada lag
    - n_pairs: número de pares para cada lag
    """
    # Calcular matriz de distâncias
    dist_matrix = squareform(pdist(coords))
    n = len(values)
    
    # Se max_dist não for especificado, usar metade da distância máxima
    if max_dist is None:
        max_dist = np.max(dist_matrix) / 2
    
    # Calcular tamanho do passo (lag)
    lag_size = max_dist / n_lags
    
    # Converter ângulos para radianos
    angle_rad = np.radians(angle)
    angle_tolerance_rad = np.radians(angle_tolerance)
    
    # Inicializar arrays para armazenar resultados
    lag_distances = np.zeros(n_lags)
    gamma = np.zeros(n_lags)
    n_pairs = np.zeros(n_lags, dtype=int)
    
    # Para cada par de pontos
    for i in range(n):
        for j in range(i+1, n):
            # Calcular distância
            dist = dist_matrix[i, j]
            
            # Verificar se está dentro da distância máxima
            if dist <= max_dist:
                # Calcular o ângulo entre os pontos
                dx = coords[j, 0] - coords[i, 0]
                dy = coords[j, 1] - coords[i, 1]
                point_angle = np.arctan2(dy, dx) % (2 * np.pi)
                
                # Verificar se o ângulo está dentro da tolerância
                angle_diff = np.abs((point_angle - angle_rad) % (2 * np.pi))
                angle_diff = min(angle_diff, 2 * np.pi - angle_diff)
                
                if angle_diff <= angle_tolerance_rad or angle_tolerance == 180:
                    # Determinar o lag
                    lag = int(dist / lag_size)
                    if lag < n_lags:
                        # Calcular a diferença quadrática
                        diff_squared = (values[i] - values[j]) ** 2
                        
                        # Atualizar os arrays
                        lag_distances[lag] += dist
                        gamma[lag] += diff_squared / 2
                        n_pairs[lag] += 1
    
    # Calcular médias para cada lag
    for lag in range(n_lags):
        if n_pairs[lag] > 0:
            lag_distances[lag] /= n_pairs[lag]
            gamma[lag] /= n_pairs[lag]
        else:
            lag_distances[lag] = np.nan
            gamma[lag] = np.nan
    
    # Remover lags sem pares
    valid_lags = ~np.isnan(lag_distances)
    lag_distances = lag_distances[valid_lags]
    gamma = gamma[valid_lags]
    n_pairs = n_pairs[valid_lags]
    
    return lag_distances, gamma, n_pairs

# Modelos teóricos de variograma
def modelo_esferico(h, c0, c1, a):
    """
    Modelo esférico de variograma.
    
    Parâmetros:
    - h: distância
    - c0: efeito pepita (nugget)
    - c1: patamar parcial (partial sill)
    - a: alcance (range)
    """
    result = np.zeros_like(h)
    mask = h <= a
    result[mask] = c0 + c1 * (1.5 * (h[mask] / a) - 0.5 * (h[mask] / a) ** 3)
    result[~mask] = c0 + c1
    return result

def modelo_exponencial(h, c0, c1, a):
    """
    Modelo exponencial de variograma.
    
    Parâmetros:
    - h: distância
    - c0: efeito pepita (nugget)
    - c1: patamar parcial (partial sill)
    - a: alcance (range)
    """
    return c0 + c1 * (1 - np.exp(-3 * h / a))

def modelo_gaussiano(h, c0, c1, a):
    """
    Modelo gaussiano de variograma.
    
    Parâmetros:
    - h: distância
    - c0: efeito pepita (nugget)
    - c1: patamar parcial (partial sill)
    - a: alcance (range)
    """
    return c0 + c1 * (1 - np.exp(-3 * h**2 / a**2))

# Função para ajustar modelos teóricos ao variograma experimental
def ajustar_modelos(lag_distances, gamma):
    """
    Ajusta os três modelos teóricos ao variograma experimental.
    
    Parâmetros:
    - lag_distances: distâncias médias para cada lag
    - gamma: valores do semivariograma para cada lag
    
    Retorna:
    - params_sph: parâmetros do modelo esférico [c0, c1, a]
    - params_exp: parâmetros do modelo exponencial [c0, c1, a]
    - params_gau: parâmetros do modelo gaussiano [c0, c1, a]
    """
    # Estimativas iniciais
    c0_init = np.min(gamma) * 0.1  # Efeito pepita inicial (10% do valor mínimo)
    c1_init = np.max(gamma) - c0_init  # Patamar parcial inicial
    a_init = np.max(lag_distances) * 0.5  # Alcance inicial (50% da distância máxima)
    
    # Limites para os parâmetros
    bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
    
    try:
        # Ajustar modelo esférico
        params_sph, _ = curve_fit(modelo_esferico, lag_distances, gamma, 
                                 p0=[c0_init, c1_init, a_init], bounds=bounds)
    except:
        params_sph = [c0_init, c1_init, a_init]
        print("Aviso: Não foi possível ajustar o modelo esférico. Usando estimativas iniciais.")
    
    try:
        # Ajustar modelo exponencial
        params_exp, _ = curve_fit(modelo_exponencial, lag_distances, gamma, 
                                 p0=[c0_init, c1_init, a_init], bounds=bounds)
    except:
        params_exp = [c0_init, c1_init, a_init]
        print("Aviso: Não foi possível ajustar o modelo exponencial. Usando estimativas iniciais.")
    
    try:
        # Ajustar modelo gaussiano
        params_gau, _ = curve_fit(modelo_gaussiano, lag_distances, gamma, 
                                 p0=[c0_init, c1_init, a_init], bounds=bounds)
    except:
        params_gau = [c0_init, c1_init, a_init]
        print("Aviso: Não foi possível ajustar o modelo gaussiano. Usando estimativas iniciais.")
    
    return params_sph, params_exp, params_gau

# Função para calcular o erro quadrático médio (MSE)
def calcular_mse(y_true, y_pred):
    """
    Calcula o erro quadrático médio entre os valores observados e preditos.
    """
    return np.mean((y_true - y_pred) ** 2)

# Definir parâmetros para a análise variográfica
print("\n=== ANÁLISE VARIOGRÁFICA ===\n")

# Parâmetros do variograma
n_lags = 10  # Número de passos (lags)
max_dist = None  # Distância máxima (será calculada como metade da distância máxima)
tolerance = 0.5  # Tolerância do passo (lag tolerance)

# Definir direções para análise (em graus, 0 = leste, 90 = norte)
direcoes = [0, 45, 90, 135]
tolerance_angular = 45  # Tolerância angular (em graus)

# Armazenar resultados para cada direção
resultados = {}

# Para cada direção
for angulo in direcoes:
    print(f"\nAnalisando direção: {angulo}°")
    
    # Calcular variograma experimental
    lag_distances, gamma, n_pairs = calcular_variograma_experimental(
        coords, values, n_lags=n_lags, max_dist=max_dist, 
        tolerance=tolerance, angle=angulo, angle_tolerance=tolerance_angular
    )
    
    # Ajustar modelos teóricos
    params_sph, params_exp, params_gau = ajustar_modelos(lag_distances, gamma)
    
    # Calcular valores preditos pelos modelos
    gamma_sph = modelo_esferico(lag_distances, *params_sph)
    gamma_exp = modelo_exponencial(lag_distances, *params_exp)
    gamma_gau = modelo_gaussiano(lag_distances, *params_gau)
    
    # Calcular MSE para cada modelo
    mse_sph = calcular_mse(gamma, gamma_sph)
    mse_exp = calcular_mse(gamma, gamma_exp)
    mse_gau = calcular_mse(gamma, gamma_gau)
    
    # Exibir parâmetros e MSE
    print(f"Modelo Esférico: C0={params_sph[0]:.4f}, C1={params_sph[1]:.4f}, a={params_sph[2]:.4f}, MSE={mse_sph:.6f}")
    print(f"Modelo Exponencial: C0={params_exp[0]:.4f}, C1={params_exp[1]:.4f}, a={params_exp[2]:.4f}, MSE={mse_exp:.6f}")
    print(f"Modelo Gaussiano: C0={params_gau[0]:.4f}, C1={params_gau[1]:.4f}, a={params_gau[2]:.4f}, MSE={mse_gau:.6f}")
    
    # Determinar o melhor modelo
    mse_values = [mse_sph, mse_exp, mse_gau]
    best_model_idx = np.argmin(mse_values)
    best_models = ['Esférico', 'Exponencial', 'Gaussiano']
    print(f"Melhor modelo: {best_models[best_model_idx]} (MSE={mse_values[best_model_idx]:.6f})")
    
    # Plotar variograma experimental e modelos ajustados
    plt.figure(figsize=(10, 6))
    
    # Plotar pontos experimentais
    plt.scatter(lag_distances, gamma, color='black', label='Experimental')
    
    # Plotar número de pares como tamanho dos pontos
    for i, (x, y, n) in enumerate(zip(lag_distances, gamma, n_pairs)):
        plt.annotate(f"{n}", (x, y), xytext=(0, 10), textcoords='offset points', ha='center')
    
    # Criar array de distâncias para plotar curvas suaves
    h = np.linspace(0, np.max(lag_distances) * 1.1, 100)
    
    # Plotar modelos ajustados
    plt.plot(h, modelo_esferico(h, *params_sph), 'r-', label=f'Esférico (MSE={mse_sph:.6f})')
    plt.plot(h, modelo_exponencial(h, *params_exp), 'g-', label=f'Exponencial (MSE={mse_exp:.6f})')
    plt.plot(h, modelo_gaussiano(h, *params_gau), 'b-', label=f'Gaussiano (MSE={mse_gau:.6f})')
    
    # Adicionar rótulos e título
    plt.title(f'Variograma Direcional - {angulo}°')
    plt.xlabel('Distância (h)')
    plt.ylabel('Semivariância γ(h)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Salvar a figura
    plt.savefig(f'figuras/variograma_{angulo}graus.png', dpi=300, bbox_inches='tight')
    
    # Armazenar resultados
    resultados[angulo] = {
        'lag_distances': lag_distances,
        'gamma': gamma,
        'n_pairs': n_pairs,
        'params_sph': params_sph,
        'params_exp': params_exp,
        'params_gau': params_gau,
        'mse_sph': mse_sph,
        'mse_exp': mse_exp,
        'mse_gau': mse_gau,
        'best_model': best_models[best_model_idx],
        'best_model_params': [params_sph, params_exp, params_gau][best_model_idx]
    }

# Determinar a melhor direção com base no menor MSE do melhor modelo para cada direção
melhor_mse = float('inf')
melhor_direcao = None
melhor_modelo = None

for angulo, res in resultados.items():
    mse_values = [res['mse_sph'], res['mse_exp'], res['mse_gau']]
    min_mse = min(mse_values)
    
    if min_mse < melhor_mse:
        melhor_mse = min_mse
        melhor_direcao = angulo
        melhor_modelo = res['best_model']

print(f"\n=== RESULTADO FINAL ===\n")
print(f"Melhor direção: {melhor_direcao}°")
print(f"Melhor modelo: {melhor_modelo}")
print(f"MSE: {melhor_mse:.6f}")

# Salvar os resultados em um arquivo
with open('resultados_variografia.txt', 'w') as f:
    f.write("=== RESULTADOS DA ANÁLISE VARIOGRÁFICA ===\n\n")
    
    for angulo, res in resultados.items():
        f.write(f"Direção: {angulo}°\n")
        f.write(f"Modelo Esférico: C0={res['params_sph'][0]:.4f}, C1={res['params_sph'][1]:.4f}, a={res['params_sph'][2]:.4f}, MSE={res['mse_sph']:.6f}\n")
        f.write(f"Modelo Exponencial: C0={res['params_exp'][0]:.4f}, C1={res['params_exp'][1]:.4f}, a={res['params_exp'][2]:.4f}, MSE={res['mse_exp']:.6f}\n")
        f.write(f"Modelo Gaussiano: C0={res['params_gau'][0]:.4f}, C1={res['params_gau'][1]:.4f}, a={res['params_gau'][2]:.4f}, MSE={res['mse_gau']:.6f}\n")
        f.write(f"Melhor modelo: {res['best_model']} (MSE={min([res['mse_sph'], res['mse_exp'], res['mse_gau']]):.6f})\n\n")
    
    f.write("=== RESULTADO FINAL ===\n")
    f.write(f"Melhor direção: {melhor_direcao}°\n")
    f.write(f"Melhor modelo: {melhor_modelo}\n")
    f.write(f"MSE: {melhor_mse:.6f}\n")

print("\nResultados salvos em 'resultados_variografia.txt'")

# Mostrar os gráficos
plt.show()
