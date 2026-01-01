FROM continuumio/miniconda3:23.5.2-0

# System deps
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy environment
COPY docker/environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -afy

# Activate env
ENV PATH /opt/conda/envs/vqc-a1/bin:$PATH

# Copy source
WORKDIR /app
COPY . /app

# Install package
RUN pip install -e .

# Default command
CMD ["vqc-drug-a1", "--help"]
