import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
from sklearn.preprocessing import PowerTransformer, QuantileTransformer

# Configurações para os gráficos
plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Carregar os dados limpos
df = pd.read_csv('dados_potassio_limpos.csv')

# Criar diretório para salvar as figuras se não existir
if not os.path.exists('figuras'):
    os.makedirs('figuras')

# 1. Análise estatística descritiva
print("\n=== ANÁLISE ESTATÍSTICA DESCRITIVA ===\n")
print("Estatísticas básicas:")
print(df['K_mg'].describe())

print("\nMedidas de assimetria e curtose:")
print(f"Assimetria: {stats.skew(df['K_mg'])}")
print(f"Curtose: {stats.kurtosis(df['K_mg'])}")

# 2. Histograma da variável
plt.figure(figsize=(10, 6))
sns.histplot(df['K_mg'], kde=True, bins=15)
plt.title('Histograma do Teor de Potássio (mg)')
plt.xlabel('Teor de Potássio (mg)')
plt.ylabel('Frequência')
plt.grid(True, alpha=0.3)
plt.savefig('figuras/histograma_potassio.png', dpi=300, bbox_inches='tight')

# Verificar assimetria e transformar os dados se necessário
assimetria = stats.skew(df['K_mg'])
print(f"\nAssimetria dos dados originais: {assimetria}")

# Se a assimetria for positiva (> 0), aplicar transformação logarítmica
if assimetria > 0:
    print("\nDetectada assimetria positiva. Aplicando transformações:")
    
    # Transformação logarítmica (adicionando 1 para evitar log(0))
    df['K_mg_log'] = np.log1p(df['K_mg'])
    
    # Transformação Box-Cox
    pt = PowerTransformer(method='box-cox')
    df['K_mg_boxcox'] = pt.fit_transform(df[['K_mg']])
    
    # Transformação Quantílica (para distribuição normal)
    qt = QuantileTransformer(output_distribution='normal')
    df['K_mg_quantile'] = qt.fit_transform(df[['K_mg']])
    
    # Plotar histogramas das transformações
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Original
    sns.histplot(df['K_mg'], kde=True, ax=axes[0, 0], bins=15)
    axes[0, 0].set_title(f'Original (Assimetria: {stats.skew(df["K_mg"]):.3f})')
    
    # Log
    sns.histplot(df['K_mg_log'], kde=True, ax=axes[0, 1], bins=15)
    axes[0, 1].set_title(f'Logarítmica (Assimetria: {stats.skew(df["K_mg_log"]):.3f})')
    
    # Box-Cox
    sns.histplot(df['K_mg_boxcox'], kde=True, ax=axes[1, 0], bins=15)
    axes[1, 0].set_title(f'Box-Cox (Assimetria: {stats.skew(df["K_mg_boxcox"]):.3f})')
    
    # Quantílica
    sns.histplot(df['K_mg_quantile'], kde=True, ax=axes[1, 1], bins=15)
    axes[1, 1].set_title(f'Quantílica (Assimetria: {stats.skew(df["K_mg_quantile"]):.3f})')
    
    plt.tight_layout()
    plt.savefig('figuras/transformacoes_potassio.png', dpi=300, bbox_inches='tight')
    
    # Salvar os dados transformados
    df.to_csv('dados_potassio_transformados.csv', index=False)
    print("\nDados transformados salvos em 'dados_potassio_transformados.csv'")

# Mostrar os gráficos
plt.show()
