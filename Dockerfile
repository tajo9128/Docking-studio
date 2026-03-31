# BioDockify Dockerfile - Full version with RDKit, OpenMM, AutoDock Vina
FROM python:3.11-slim

LABEL maintainer="BioDockify"
LABEL description="Molecular Docking Studio with AI Assistant"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update && apt-get install -y \
    wget \
    curl \
    gnupg2 \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    fastapi==0.109.2 \
    uvicorn[standard]==0.27.1 \
    python-multipart==0.0.9 \
    pydantic==2.6.1 \
    pydantic-settings==2.1.0 \
    httpx==0.26.0 \
    numpy==1.26.4 \
    pandas==2.2.1

RUN pip install --no-cache-dir rdkit openmm

RUN apt-get update && apt-get install -y autodock-vina openbabel && rm -rf /var/lib/apt/lists/*

COPY app.py /app/

RUN mkdir -p /app/data/jobs

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8000/health || exit 1

CMD ["python", "/app/app.py"]
