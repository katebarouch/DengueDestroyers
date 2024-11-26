#!/bin/bash

### REQUIREMENTS ###
# Use this file to set up all necessary dependencies for running the project on Linux (Ubuntu/Debian).

### Step 1: Install Necessary Tools ###

# 1.1 Install wget
sudo apt install wget

# 1.2 Install samtools
sudo apt install samtools

# 1.3 Install bowtie2
sudo apt install bowtie2

# 1.4 Install Node.js and JBrowse CLI
sudo apt install nodejs
sudo npm install -g @jbrowse/cli

# Verify installation:
node -v
jbrowse --version

### Step 2: Apache Server Setup ###

# 2.1 Install apache2
sudo apt install apache2

# 2.2 Start the apache2 server
sudo service apache2 start

### Additional Resources ###
# - Apache2 setup: https://httpd.apache.org/docs/
# - JBrowse CLI: https://jbrowse.org/
