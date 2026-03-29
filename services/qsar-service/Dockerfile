FROM python:3.11-slim

LABEL maintainer="BioDockify"
LABEL description="QSAR Service - ML-based QSAR modeling"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    libopenbabel7 \
    libhdf5-dev \
    libopenblas-dev \
    libglib2.0-0 \
    libgfortran5 \
    libsm6 \
    libxml2 \
    libxslt1.1 \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    rdkit-pypi \
    fastapi \
    uvicorn \
    pydantic \
    python-multipart==0.0.22 \
    requests \
    redis \
    numpy \
    pandas \
    scikit-learn \
    plotly \
    matplotlib \
    joblib \
    httpx

WORKDIR /app

COPY . .

EXPOSE 8005

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8005"]
