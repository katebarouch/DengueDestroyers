#!/bin/bash

### Dengue Alignments Tracks ###

# Helper function for error handling
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error encountered in previous command. Exiting."
        exit 1
    fi
}

### Check Directory and APACHE_ROOT ###
WORKDIR=$(pwd)

# Check if APACHE_ROOT is defined
if [ -z "$APACHE_ROOT" ]; then
    echo "Error: APACHE_ROOT is not defined. Please set it using the following command:"
    echo "export APACHE_ROOT='/path/to/apache/root'"
    echo "Refer to the instructions in the script header or requirements.txt for details."
    exit 1
fi

echo "Working directory set to: $WORKDIR"
echo "Using Apache root directory: $APACHE_ROOT"

# Ensure JBrowse2 directory exists
if [ ! -d "$APACHE_ROOT/jbrowse2" ]; then
    echo "Error: JBrowse2 not found in $APACHE_ROOT/jbrowse2."
    echo "Please ensure JBrowse2 is installed and configured correctly."
    exit 1
fi

### Download and Process the Reference Genomes for all serotypes ###

# Define URLs for reference genomes
declare -A genome_urls # -A is for Associative (aka string array)
genome_urls["DENV-1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz"
genome_urls["DENV-2"]="..."
genome_urls["DENV-3"]="..."
genome_urls["DENV-4"]="..."

# Define URLs for reference annotations
declare -A annotation_urls
annotation_urls["DENV-1"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz"
annotation_urls["DENV-2"]="..."
annotation_urls["DENV-3"]="..."
annotation_urls["DENV-4"]="..."

# Define URLs and JBrowse names for comparison genomes (all 2 entries for each serotype are stored in single string for simplicity)
declare -A synteny_comparison_genomes
comparison_genomes["DENV-1"]="https://example.com/DENV_1_comparison1.fa,SampleA https://example.com/DENV_1_comparison2.fa,SampleB"
comparison_genomes["DENV-2"]="https://example.com/DENV_2_comparison1.fa,IsolateX https://example.com/DENV_2_comparison2.fa,IsolateY"
comparison_genomes["DENV-3"]="https://example.com/DENV_3_comparison1.fa,Variant1 https://example.com/DENV_3_comparison2.fa,Variant2"
comparison_genomes["DENV-4"]="https://example.com/DENV_4_comparison1.fa,StrainAlpha https://example.com/DENV_4_comparison2.fa,StrainBeta"

