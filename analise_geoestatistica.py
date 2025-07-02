import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Configurações para os gráficos
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Função para carregar os dados
def carregar_dados(arquivo):
    # Ler o arquivo manualmente para lidar com o formato específico
    with open(arquivo, 'r') as f:
        linhas = f.readlines()
    
    # Pular as primeiras 5 linhas (cabeçalho)
    dados = []
    for linha in linhas[5:]:  # Começar da linha 6 (índice 5)
        valores = linha.strip().split()
        if len(valores) == 3:
            try:
                x, y, k = valores
                # Converter para float, mas tratar valores especiais
                x_val = float(x)
                y_val = float(y)
                k_val = float(k) if k != '-999' else np.nan
                dados.append([x_val, y_val, k_val])
            except ValueError:
                # Pular linhas que não podem ser convertidas para float
                continue
    
    df = pd.DataFrame(dados, columns=['x', 'y', 'K_mg'])
    
    # Remover valores NaN
    df = df.dropna()
    
    return df

# Carregar os dados
arquivo = 'Teor_potassio (mg).txt'
df = carregar_dados(arquivo)

# Exibir informações básicas sobre os dados
print("Informações sobre o conjunto de dados:")
print(f"Número de amostras: {len(df)}")
print("\nEstatísticas descritivas:")
print(df.describe())

# Salvar o DataFrame limpo para uso posterior
df.to_csv('dados_potassio_limpos.csv', index=False)

print("\nDados limpos salvos em 'dados_potassio_limpos.csv'")
