#!/bin/bash

### REQUIREMENTS ###
# Use this file to set up all necessary dependencies for running the project on macOS.

### Step 1: Install Necessary Tools ###

# 1.1 Install wget
brew install wget

# 1.2 Install samtools
brew install samtools

# 1.3 Install bowtie2
brew install bowtie2

# 1.4 Install Node.js and JBrowse CLI
brew install node@20
sudo npm install -g @jbrowse/cli

# Verify installation:
node -v
jbrowse --version

### Step 2: Apache Server Setup ###

# 2.1 Install apache2
brew install httpd

# 2.2 Start the apache2 server
sudo brew services start httpd

### Additional Resources ###
# - Apache2 setup: https://httpd.apache.org/docs/
# - JBrowse CLI: https://jbrowse.org/
