FROM python:3.14.7-slim-bookworm

# Prevent Python from writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1

# Prevent Python output from being buffered
ENV PYTHONUNBUFFERED=1

WORKDIR /app

# System dependencies required by scientific Python packages
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    libhdf5-dev \
    libgl1 \
    libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies first for Docker layer caching
COPY requirements.txt .

RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copy application
COPY . .

# Streamlit configuration
ENV STREAMLIT_SERVER_HEADLESS=true
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.maxUploadSize=500"]