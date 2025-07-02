import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.spatial.distance import pdist, squareform
import os
from sklearn.preprocessing import PowerTransformer
from matplotlib.colors import Normalize
from matplotlib import cm

# Configurações para os gráficos
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Carregar os dados limpos
df = pd.read_csv('dados_potassio_limpos.csv')

# Criar diretório para salvar as figuras se não existir
if not os.path.exists('figuras'):
    os.makedirs('figuras')

# 1. Visualização da distribuição espacial dos pontos
plt.figure(figsize=(10, 8))

# Criar um scatter plot com cores baseadas no valor de K_mg
scatter = plt.scatter(df['x'], df['y'], c=df['K_mg'], cmap='viridis', 
                     s=100, alpha=0.8, edgecolors='k', linewidths=0.5)

# Adicionar barra de cores
cbar = plt.colorbar(scatter)
cbar.set_label('Teor de Potássio (mg)')

# Adicionar rótulos e título
plt.title('Distribuição Espacial do Teor de Potássio')
plt.xlabel('Coordenada X')
plt.ylabel('Coordenada Y')
plt.grid(True, alpha=0.3)

# Salvar a figura
plt.savefig('figuras/distribuicao_espacial.png', dpi=300, bbox_inches='tight')

# 2. Definir parâmetros para a grid
print("\n=== DEFINIÇÃO DA GRID ===\n")

# Calcular os limites dos dados
x_min, x_max = df['x'].min(), df['x'].max()
y_min, y_max = df['y'].min(), df['y'].max()

# Adicionar uma margem de 5% aos limites
margem_x = 0.05 * (x_max - x_min)
margem_y = 0.05 * (y_max - y_min)

x_min -= margem_x
x_max += margem_x
y_min -= margem_y
y_max += margem_y

# Definir o tamanho das células da grid (10x10 células)
tamanho_celula_x = (x_max - x_min) / 10
tamanho_celula_y = (y_max - y_min) / 10

# Arredondar para facilitar a interpretação
tamanho_celula_x = round(tamanho_celula_x, 1)
tamanho_celula_y = round(tamanho_celula_y, 1)

# Recalcular os limites para ter células de tamanho exato
x_max = x_min + 10 * tamanho_celula_x
y_max = y_min + 10 * tamanho_celula_y

# Criar a grid
x_grid = np.arange(x_min, x_max + tamanho_celula_x, tamanho_celula_x)
y_grid = np.arange(y_min, y_max + tamanho_celula_y, tamanho_celula_y)

# Exibir informações da grid
print(f"Dimensão da grid: 10 x 10 células")
print(f"Tamanho das células: {tamanho_celula_x} x {tamanho_celula_y} unidades")
print(f"Coordenada de origem: ({x_min}, {y_min})")
print(f"Coordenada final: ({x_max}, {y_max})")

# 3. Visualizar a grid com os pontos
plt.figure(figsize=(10, 8))

# Plotar a grid
for x in x_grid:
    plt.axvline(x=x, color='gray', linestyle='-', alpha=0.3)
for y in y_grid:
    plt.axhline(y=y, color='gray', linestyle='-', alpha=0.3)

# Plotar os pontos
scatter = plt.scatter(df['x'], df['y'], c=df['K_mg'], cmap='viridis', 
                     s=100, alpha=0.8, edgecolors='k', linewidths=0.5)

# Adicionar barra de cores
cbar = plt.colorbar(scatter)
cbar.set_label('Teor de Potássio (mg)')

# Adicionar rótulos e título
plt.title('Grid Ajustada com Pontos Amostrais')
plt.xlabel('Coordenada X')
plt.ylabel('Coordenada Y')
plt.xlim(x_min - margem_x, x_max + margem_x)
plt.ylim(y_min - margem_y, y_max + margem_y)

# Salvar a figura
plt.savefig('figuras/grid_com_pontos.png', dpi=300, bbox_inches='tight')

# 4. Preparar para análise variográfica
print("\n=== PREPARAÇÃO PARA ANÁLISE VARIOGRÁFICA ===\n")

# Extrair coordenadas e valores
coords = df[['x', 'y']].values
values = df['K_mg'].values

# Salvar os parâmetros da grid para uso posterior
parametros_grid = {
    'x_min': x_min,
    'x_max': x_max,
    'y_min': y_min,
    'y_max': y_max,
    'tamanho_celula_x': tamanho_celula_x,
    'tamanho_celula_y': tamanho_celula_y,
    'n_celulas_x': 10,
    'n_celulas_y': 10
}

# Salvar os parâmetros em um arquivo
with open('parametros_grid.txt', 'w') as f:
    for key, value in parametros_grid.items():
        f.write(f"{key}: {value}\n")

print("Parâmetros da grid salvos em 'parametros_grid.txt'")

# Mostrar os gráficos
plt.show()
