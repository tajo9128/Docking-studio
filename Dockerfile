# BioDockify Docking Studio - All-in-One Docker Image
# Auto-detects GPU and uses appropriate engine (GNINA for GPU, Vina for CPU)

# Stage 1: GNINA build (from pre-built image)
FROM gnina/gnina:latest as gnina-stage

# Stage 2: Final Image (CPU + GPU)
FROM python:3.9-slim

LABEL maintainer="BioDockify"
LABEL description="BioDockify Docking Studio - Auto-detects GPU and uses GNINA or Vina"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    libhdf5-serial-103 \
    libopenblas-base \
    libglib2.0-0 \
    libgfortran5 \
    libsm6 \
    libxml2 \
    libxslt1.1 \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install Vina via pip
RUN pip install --no-cache-dir vina

# Copy GNINA from gnina stage
COPY --from=gnina-stage /usr/local/bin/gnina /usr/local/bin/gnina

# Set working directory
WORKDIR /data

# Copy entrypoint
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENV RECEPTOR_FILE=/data/receptor.pdbqt
ENV LIGAND_FILE=/data/ligand.pdbqt

ENTRYPOINT ["/entrypoint.sh"]
