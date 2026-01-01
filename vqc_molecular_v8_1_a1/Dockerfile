FROM continuumio/miniconda3:23.5.2-0

# Metadados
LABEL maintainer="Quantum Drug Discovery Team"
LABEL description="VQC-Molecular v8.1-A1: Quantum-Enhanced Drug Screening with Qualis A1 Statistical Auditing"
LABEL version="8.1-A1"

# Copiar environment.yml
COPY environment.yml .

# Criar environment
RUN conda env create -f environment.yml && \
    conda clean -afy && \
    echo "source activate vqc-a1" > ~/.bashrc

# Definir PATH
ENV PATH /opt/conda/envs/vqc-a1/bin:$PATH
ENV PYTHONUNBUFFERED=1

# Copiar código
COPY . /app
WORKDIR /app

# Criar diretórios
RUN mkdir -p results_vqc_drug qsar_cache logs 01_protocolo_pre_registrado \
    02_dados_brutos 03_analises_estatisticas 04_figuras_publicacao \
    05_tabelas_suplementares 06_reprodutibilidade 07_log_execucao

# Comando padrão
CMD ["python", "vqc_drug_qualis_a1.py", "--target", "EGFR", "--trials", "500"]
