#!/bin/bash

### REQUIREMENTS ###
# Use this file to set up all necessary dependencies for running the project on Linux (Ubuntu/Debian).

### Step 1: Install Necessary Tools ###

# 1.1 Install wget
if ! command -v wget &> /dev/null; then
    echo "wget is not installed. Installing now..."
    sudo apt update && sudo apt install -y wget
else
    echo "wget is already installed."
fi

# 1.2 Install samtools
if ! command -v samtools &> /dev/null; then
    echo "samtools is not installed. Installing now..."
    sudo apt install -y samtools
else
    echo "samtools is already installed."
fi

# 1.3 Install bowtie2
if ! command -v bowtie2 &> /dev/null; then
    echo "bowtie2 is not installed. Installing now..."
    sudo apt install -y bowtie2
else
    echo "bowtie2 is already installed."
fi

# 1.4 Install Node.js
if ! command -v node &> /dev/null; then
    echo "Node.js is not installed. Installing now..."
    sudo apt install -y nodejs npm
else
    echo "Node.js is already installed."
fi

# Install JBrowse CLI
if ! command -v jbrowse &> /dev/null; then
    echo "JBrowse CLI is not installed. Installing now..."
    sudo npm install -g @jbrowse/cli
else
    echo "JBrowse CLI is already installed."
fi
# Verify installation:
node -v
jbrowse --version

### Step 2: Apache Server Setup ###

# 2.1 Install apache2
echo "Installing apache2..."
sudo apt install apache2

# 2.2 Start the apache2 server
echo "Starting apache2 server..."
sudo service apache2 start

### Additional Resources ###
# - Apache2 setup: https://httpd.apache.org/docs/
# - JBrowse CLI: https://jbrowse.org/
