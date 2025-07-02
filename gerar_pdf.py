import os
import subprocess

# Verificar se o pandoc está instalado
def verificar_pandoc():
    try:
        subprocess.run(['pandoc', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False

# Função para converter markdown para PDF
def converter_md_para_pdf(arquivo_md, arquivo_pdf):
    if not verificar_pandoc():
        print("Erro: Pandoc não está instalado. Por favor, instale o Pandoc para gerar o PDF.")
        print("Você pode instalar com: sudo apt-get install pandoc texlive-xetex")
        return False
    
    try:
        # Comando para converter markdown para PDF usando pandoc
        cmd = [
            'pandoc',
            arquivo_md,
            '-o', arquivo_pdf,
            '--pdf-engine=xelatex',
            '-V', 'geometry:margin=1in',
            '-V', 'colorlinks=true',
            '-V', 'linkcolor=blue',
            '-V', 'urlcolor=blue',
            '--toc',  # Adicionar tabela de conteúdo
            '--toc-depth=3',
            '--standalone',
            '--embed-resources',  # Incorporar recursos (imagens)
            '--resource-path=.:./figuras',  # Procurar recursos no diretório atual e na pasta figuras
            '-V', 'papersize=a4',
            '-V', 'fontsize=11pt'
        ]
        
        subprocess.run(cmd, check=True)
        print(f"PDF gerado com sucesso: {arquivo_pdf}")
        return True
    except subprocess.SubprocessError as e:
        print(f"Erro ao gerar PDF: {e}")
        return False

# Arquivo de entrada e saída
arquivo_md = 'relatorio_final.md'
arquivo_pdf = 'Relatorio_Final_Geoestatistica_Pedro_Henrique_Silva_18101395.pdf'

# Verificar se o arquivo markdown existe
if not os.path.exists(arquivo_md):
    print(f"Erro: O arquivo {arquivo_md} não foi encontrado.")
else:
    # Converter para PDF
    converter_md_para_pdf(arquivo_md, arquivo_pdf)
