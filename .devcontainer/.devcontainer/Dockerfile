# Minimal devcontainer with BLAST, Python3, pip, and R
FROM mcr.microsoft.com/devcontainers/base:ubuntu

# Install system packages, BLAST, R, Python, and common build headers
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    r-base \
    python3 \
    python3-pip \
    build-essential \
    wget \
    curl \
    git \
    tree \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Make python and pip point to python3/pip3
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1 \
 && update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

WORKDIR /workspace
