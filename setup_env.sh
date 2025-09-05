#!/bin/bash

# Update system
sudo apt-get update && sudo apt-get install -y \
    python3-pip \
    r-base \
    ncbi-blast+ \
    build-essential \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    wget \
    curl \
    git

# Python packages
pip3 install --user biopython pandas numpy matplotlib seaborn

# Optional: R packages
Rscript -e "install.packages(c('tidyverse'), repos='https://cloud.r-project.org')"

echo "Environment setup complete!"
