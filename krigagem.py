import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import os
from scipy.spatial.distance import cdist

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

# Carregar parâmetros da grid
parametros = {}
with open('parametros_grid.txt', 'r') as f:
    for linha in f:
        chave, valor = linha.strip().split(': ')
        parametros[chave] = float(valor)

# Definir a grid para estimativa
x_min = parametros['x_min']
x_max = parametros['x_max']
y_min = parametros['y_min']
y_max = parametros['y_max']
tamanho_celula_x = parametros['tamanho_celula_x']
tamanho_celula_y = parametros['tamanho_celula_y']

# Criar grid de pontos para estimativa
x_grid = np.arange(x_min, x_max + tamanho_celula_x, tamanho_celula_x)
y_grid = np.arange(y_min, y_max + tamanho_celula_y, tamanho_celula_y)
X, Y = np.meshgrid(x_grid, y_grid)
grid_points = np.column_stack([X.flatten(), Y.flatten()])

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

# Função para calcular a matriz de covariância
def calcular_matriz_covariancia(coords1, coords2, modelo, params):
    """
    Calcula a matriz de covariância entre dois conjuntos de coordenadas.
    
    Parâmetros:
    - coords1: primeiro conjunto de coordenadas
    - coords2: segundo conjunto de coordenadas
    - modelo: função do modelo de variograma
    - params: parâmetros do modelo [c0, c1, a]
    
    Retorna:
    - matriz de covariância
    """
    # Calcular matriz de distâncias
    dist_matrix = cdist(coords1, coords2)
    
    # Calcular matriz de variograma
    variogram_matrix = modelo(dist_matrix, *params)
    
    # Converter para matriz de covariância
    c0, c1, _ = params
    sill = c0 + c1  # Patamar total
    cov_matrix = sill - variogram_matrix
    
    return cov_matrix

# Função para realizar krigagem ordinária
def krigagem_ordinaria(coords_amostra, valores_amostra, coords_estimativa, modelo, params):
    """
    Realiza krigagem ordinária.
    
    Parâmetros:
    - coords_amostra: coordenadas dos pontos amostrais
    - valores_amostra: valores dos pontos amostrais
    - coords_estimativa: coordenadas dos pontos a estimar
    - modelo: função do modelo de variograma
    - params: parâmetros do modelo [c0, c1, a]
    
    Retorna:
    - estimativas: valores estimados
    - variancias: variâncias de krigagem
    """
    n_amostras = len(coords_amostra)
    n_estimativas = len(coords_estimativa)
    
    # Inicializar arrays para armazenar resultados
    estimativas = np.zeros(n_estimativas)
    variancias = np.zeros(n_estimativas)
    
    # Calcular matriz de covariância entre amostras
    K = calcular_matriz_covariancia(coords_amostra, coords_amostra, modelo, params)
    
    # Adicionar restrição de soma dos pesos = 1 (multiplicador de Lagrange)
    K_extended = np.zeros((n_amostras + 1, n_amostras + 1))
    K_extended[:n_amostras, :n_amostras] = K
    K_extended[n_amostras, :n_amostras] = 1
    K_extended[:n_amostras, n_amostras] = 1
    
    # Para cada ponto a estimar
    for i in range(n_estimativas):
        # Calcular covariâncias entre o ponto a estimar e as amostras
        k = calcular_matriz_covariancia(np.array([coords_estimativa[i]]), coords_amostra, modelo, params)[0]
        
        # Montar vetor do lado direito
        b = np.zeros(n_amostras + 1)
        b[:n_amostras] = k
        b[n_amostras] = 1
        
        try:
            # Resolver sistema de equações para obter pesos
            weights = np.linalg.solve(K_extended, b)
            
            # Extrair pesos de krigagem (sem o multiplicador de Lagrange)
            lambda_weights = weights[:n_amostras]
            
            # Calcular estimativa
            estimativas[i] = np.sum(lambda_weights * valores_amostra)
            
            # Calcular variância de krigagem
            variancias[i] = np.sum(lambda_weights * k) + weights[n_amostras]
        except np.linalg.LinAlgError:
            # Em caso de erro, usar média simples
            estimativas[i] = np.mean(valores_amostra)
            variancias[i] = np.var(valores_amostra)
    
    return estimativas, variancias

# Carregar os parâmetros dos modelos da melhor direção (135°)
print("\n=== ESTIMATIVA POR KRIGAGEM ORDINÁRIA ===\n")

# Parâmetros dos modelos (obtidos da análise variográfica)
params_sph = [6.5572, 15.6714, 104.1889]  # [c0, c1, a] para modelo esférico
params_exp = [5.2636, 18.8732, 151.6199]  # [c0, c1, a] para modelo exponencial
params_gau = [9.1771, 13.3829, 95.2249]   # [c0, c1, a] para modelo gaussiano

print("Realizando krigagem com modelo esférico...")
# Realizar krigagem com modelo esférico
estimativas_sph, variancias_sph = krigagem_ordinaria(
    coords, values, grid_points, modelo_esferico, params_sph
)

print("Realizando krigagem com modelo exponencial...")
# Realizar krigagem com modelo exponencial
estimativas_exp, variancias_exp = krigagem_ordinaria(
    coords, values, grid_points, modelo_exponencial, params_exp
)

print("Realizando krigagem com modelo gaussiano...")
# Realizar krigagem com modelo gaussiano
estimativas_gau, variancias_gau = krigagem_ordinaria(
    coords, values, grid_points, modelo_gaussiano, params_gau
)

# Reorganizar resultados em formato de grid
shape = (len(y_grid), len(x_grid))

# Estimativas
Z_sph = estimativas_sph.reshape(shape)
Z_exp = estimativas_exp.reshape(shape)
Z_gau = estimativas_gau.reshape(shape)

# Variâncias
V_sph = variancias_sph.reshape(shape)
V_exp = variancias_exp.reshape(shape)
V_gau = variancias_gau.reshape(shape)

# Plotar mapas de estimativas
fig, axes = plt.subplots(3, 2, figsize=(18, 24))

# Definir limites de cores para estimativas e variâncias
vmin_est = min(np.min(Z_sph), np.min(Z_exp), np.min(Z_gau))
vmax_est = max(np.max(Z_sph), np.max(Z_exp), np.max(Z_gau))

vmin_var = min(np.min(V_sph), np.min(V_exp), np.min(V_gau))
vmax_var = max(np.max(V_sph), np.max(V_exp), np.max(V_gau))

# Modelo esférico
im1 = axes[0, 0].imshow(Z_sph, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                      aspect='equal', cmap='viridis', vmin=vmin_est, vmax=vmax_est)
axes[0, 0].scatter(coords[:, 0], coords[:, 1], c='red', s=30, edgecolors='k')
axes[0, 0].set_title('Estimativa - Modelo Esférico')
axes[0, 0].set_xlabel('Coordenada X')
axes[0, 0].set_ylabel('Coordenada Y')
fig.colorbar(im1, ax=axes[0, 0], label='Teor de Potássio (mg)')

im2 = axes[0, 1].imshow(V_sph, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                      aspect='equal', cmap='plasma', vmin=vmin_var, vmax=vmax_var)
axes[0, 1].scatter(coords[:, 0], coords[:, 1], c='red', s=30, edgecolors='k')
axes[0, 1].set_title('Variância - Modelo Esférico')
axes[0, 1].set_xlabel('Coordenada X')
axes[0, 1].set_ylabel('Coordenada Y')
fig.colorbar(im2, ax=axes[0, 1], label='Variância')

# Modelo exponencial
im3 = axes[1, 0].imshow(Z_exp, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                      aspect='equal', cmap='viridis', vmin=vmin_est, vmax=vmax_est)
axes[1, 0].scatter(coords[:, 0], coords[:, 1], c='red', s=30, edgecolors='k')
axes[1, 0].set_title('Estimativa - Modelo Exponencial')
axes[1, 0].set_xlabel('Coordenada X')
axes[1, 0].set_ylabel('Coordenada Y')
fig.colorbar(im3, ax=axes[1, 0], label='Teor de Potássio (mg)')

im4 = axes[1, 1].imshow(V_exp, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                      aspect='equal', cmap='plasma', vmin=vmin_var, vmax=vmax_var)
axes[1, 1].scatter(coords[:, 0], coords[:, 1], c='red', s=30, edgecolors='k')
axes[1, 1].set_title('Variância - Modelo Exponencial')
axes[1, 1].set_xlabel('Coordenada X')
axes[1, 1].set_ylabel('Coordenada Y')
fig.colorbar(im4, ax=axes[1, 1], label='Variância')

# Modelo gaussiano
im5 = axes[2, 0].imshow(Z_gau, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                      aspect='equal', cmap='viridis', vmin=vmin_est, vmax=vmax_est)
axes[2, 0].scatter(coords[:, 0], coords[:, 1], c='red', s=30, edgecolors='k')
axes[2, 0].set_title('Estimativa - Modelo Gaussiano')
axes[2, 0].set_xlabel('Coordenada X')
axes[2, 0].set_ylabel('Coordenada Y')
fig.colorbar(im5, ax=axes[2, 0], label='Teor de Potássio (mg)')

im6 = axes[2, 1].imshow(V_gau, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                      aspect='equal', cmap='plasma', vmin=vmin_var, vmax=vmax_var)
axes[2, 1].scatter(coords[:, 0], coords[:, 1], c='red', s=30, edgecolors='k')
axes[2, 1].set_title('Variância - Modelo Gaussiano')
axes[2, 1].set_xlabel('Coordenada X')
axes[2, 1].set_ylabel('Coordenada Y')
fig.colorbar(im6, ax=axes[2, 1], label='Variância')

plt.tight_layout()
plt.savefig('figuras/krigagem_comparacao.png', dpi=300, bbox_inches='tight')

# Salvar os resultados em um arquivo
with open('resultados_krigagem.txt', 'w') as f:
    f.write("=== RESULTADOS DA KRIGAGEM ORDINÁRIA ===\n\n")
    
    f.write("Parâmetros do modelo esférico:\n")
    f.write(f"C0 (efeito pepita): {params_sph[0]:.4f}\n")
    f.write(f"C1 (patamar parcial): {params_sph[1]:.4f}\n")
    f.write(f"a (alcance): {params_sph[2]:.4f}\n\n")
    
    f.write("Parâmetros do modelo exponencial:\n")
    f.write(f"C0 (efeito pepita): {params_exp[0]:.4f}\n")
    f.write(f"C1 (patamar parcial): {params_exp[1]:.4f}\n")
    f.write(f"a (alcance): {params_exp[2]:.4f}\n\n")
    
    f.write("Parâmetros do modelo gaussiano:\n")
    f.write(f"C0 (efeito pepita): {params_gau[0]:.4f}\n")
    f.write(f"C1 (patamar parcial): {params_gau[1]:.4f}\n")
    f.write(f"a (alcance): {params_gau[2]:.4f}\n\n")
    
    f.write("Estatísticas das estimativas:\n")
    f.write(f"Modelo Esférico: min={np.min(Z_sph):.4f}, max={np.max(Z_sph):.4f}, média={np.mean(Z_sph):.4f}\n")
    f.write(f"Modelo Exponencial: min={np.min(Z_exp):.4f}, max={np.max(Z_exp):.4f}, média={np.mean(Z_exp):.4f}\n")
    f.write(f"Modelo Gaussiano: min={np.min(Z_gau):.4f}, max={np.max(Z_gau):.4f}, média={np.mean(Z_gau):.4f}\n\n")
    
    f.write("Estatísticas das variâncias:\n")
    f.write(f"Modelo Esférico: min={np.min(V_sph):.4f}, max={np.max(V_sph):.4f}, média={np.mean(V_sph):.4f}\n")
    f.write(f"Modelo Exponencial: min={np.min(V_exp):.4f}, max={np.max(V_exp):.4f}, média={np.mean(V_exp):.4f}\n")
    f.write(f"Modelo Gaussiano: min={np.min(V_gau):.4f}, max={np.max(V_gau):.4f}, média={np.mean(V_gau):.4f}\n")

print("\nResultados salvos em 'resultados_krigagem.txt'")
print("Mapas de krigagem salvos em 'figuras/krigagem_comparacao.png'")

# Mostrar os gráficos
plt.show()
