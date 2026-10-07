FROM python:3.14.7-slim-bookworm

# Prevent Python from writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1

# Prevent Python output from being buffered
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# System dependencies required by scientific Python packages
RUN apt-get update && apt-get install --no-install-recommends -y \
    build-essential \
    gcc \
    g++ \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies first for Docker layer caching
COPY requirements.txt .

RUN python -c "from pathlib import Path; p = Path('requirements.txt'); p.write_text(p.read_text(encoding='utf-16'), encoding='utf-8')" \
    && pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copy only files needed to run the application
COPY app.py ./
COPY src/ ./src/
COPY pages/ ./pages/
COPY config/ ./config/
COPY docs/ ./docs/

# Run the application as an unprivileged user
RUN groupadd --system bioworkbench \
    && useradd --system --gid bioworkbench --create-home bioworkbench \
    && chown -R bioworkbench:bioworkbench /app
USER bioworkbench

# Streamlit configuration
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0
ENV STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
    CMD ["python", "-c", "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8501/_stcore/health', timeout=3)"]

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.maxUploadSize=1000"]
