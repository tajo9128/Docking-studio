# BioDockify Docking Studio - All-in-One Docker Image
# Using supervisor approach similar to Agent Zero for stable startup

# Stage 1: GNINA build (from pre-built image)
FROM gnina/gnina:latest AS gnina-stage

# Stage 2: Final Image (CPU + GPU)
FROM python:3.11-slim

LABEL maintainer="BioDockify"
LABEL description="Docking Studio - GPU auto-detection with Vina+GNINA+RF+ODDT pipeline"

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    supervisor \
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

# Install Python dependencies in correct order
RUN pip install --no-cache-dir \
    six && \
    pip install --no-cache-dir \
    vina \
    rdkit \
    meeko \
    oddt

# Copy GNINA from gnina stage
COPY --from=gnina-stage /usr/local/bin/gnina /usr/local/bin/gnina

# Set working directory
WORKDIR /app

# Copy backend files
COPY backend/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY backend/ /app/

# Create supervisor configuration
RUN mkdir -p /var/log/supervisor

# Create startup script
RUN echo '#!/bin/bash' > /startup.sh && \
    echo 'echo ""' >> /startup.sh && \
    echo 'echo "============================================================"' >> /startup.sh && \
    echo 'echo "  🧬 Docking Studio - Backend Starting..."' >> /startup.sh && \
    echo 'echo "============================================================"' >> /startup.sh && \
    echo 'echo ""' >> /startup.sh && \
    echo 'echo "  📚 API Documentation: http://localhost:8000/docs"' >> /startup.sh && \
    echo 'echo "  📖 ReDoc:            http://localhost:8000/redoc"' >> /startup.sh && \
    echo 'echo "  ✅ Health:           http://localhost:8000/health"' >> /startup.sh && \
    echo 'echo "============================================================"' >> /startup.sh && \
    echo 'echo ""' >> /startup.sh && \
    chmod +x /startup.sh

# Create supervisor config
RUN echo '[supervisord]' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'nodaemon=true' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'user=root' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'logfile=/dev/stdout' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'logfile_maxbytes=0' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'pidfile=/var/run/supervisord.pid' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo '' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo '[program:uvicorn]' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'command=uvicorn main:app --host 0.0.0.0 --port 8000 --reload' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'directory=/app' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'user=root' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'stdout_logfile=/dev/stdout' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'stdout_logfile_maxbytes=0' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'stderr_logfile=/dev/stderr' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'stderr_logfile_maxbytes=0' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'autorestart=true' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'startretries=3' >> /etc/supervisor/conf.d/docking-studio.conf && \
    echo 'stopwaitsecs=30' >> /etc/supervisor/conf.d/docking-studio.conf

# Expose port
EXPOSE 8000

# Run supervisor (which manages uvicorn)
CMD ["/bin/bash", "-c", "/startup.sh && /usr/bin/supervisord -c /etc/supervisor/conf.d/docking-studio.conf"]
