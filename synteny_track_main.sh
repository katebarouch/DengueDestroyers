#!/bin/bash

### Dengue Synteny Tracks ###

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

# Function to process a reference genome
process_reference_genome() {
    virus_name="$1"
    genome_url="$2"
    echo "Processing $virus_name reference genome..."
    wget -q "$genome_url" -O "$WORKDIR/${virus_name}_genome.fna.gz"
    check_error

    echo "Unzipping $virus_name reference genome..."
    gunzip -f "$WORKDIR/${virus_name}_genome.fna.gz"
    check_error

    echo "Indexing $virus_name reference genome..."
    samtools faidx "$WORKDIR/${virus_name}_genome.fna"
    check_error

    echo "Adding $virus_name reference genome to JBrowse..."
    jbrowse add-assembly "$WORKDIR/${virus_name}_genome.fna" --out "$APACHE_ROOT/jbrowse2" --load copy
    check_error
}

# Function to process a reference annotation
process_reference_annotation() {
    virus_name="$1"
    annotation_url="$2"
    echo "Processing $virus_name annotations..."
    wget -q "$annotation_url" -O "$WORKDIR/${virus_name}_annotations.gff.gz"
    check_error

    echo "Unzipping $virus_name annotations..."
    gunzip -f "$WORKDIR/${virus_name}_annotations.gff.gz"
    check_error

    echo "Sorting $virus_name annotations..."
    jbrowse sort-gff "$WORKDIR/${virus_name}_annotations.gff" > "$WORKDIR/${virus_name}_genes.gff"
    check_error

    echo "Compressing $virus_name annotations..."
    bgzip -f "$WORKDIR/${virus_name}_genes.gff"
    check_error

    echo "Indexing $virus_name annotations..."
    tabix "$WORKDIR/${virus_name}_genes.gff.gz"
    check_error

    echo "Adding $virus_name annotations to JBrowse..."
    jbrowse add-track "$WORKDIR/${virus_name}_genes.gff.gz" --out "$APACHE_ROOT/jbrowse2" --load copy
    check_error
}

# # Function to process comparison genome for synteny track
# process_comparison_genome() {
#     reference_virus_name="$1"
#     comparison_virus_name="$2"
#     comparison_genome_url="$3"

#    ...
# }

### Main Script Execution ###

# Process DENV-1 Genomes
process_reference_genome "DENV-1" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz"
process_reference_annotation "DENV-1" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz"

process_comparison_genome "DENV-1" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" \
    "Variant1" "https://example.com/DENV_1_comparison1.fa"
process_comparison_genome "DENV-1" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" \
    "Variant2" "https://example.com/DENV_1_comparison2.fa"
process_comparison_genome "DENV-1" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" \
    "Variant3" "https://example.com/DENV_1_comparison3.fa"

# Process DENV-2 Genomes
process_reference_genome "DENV-2" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/871/845/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.fna.gz"
process_reference_annotation "DENV-2" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/871/845/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.gff.gz"

process_comparison_genome "DENV-2" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/871/845/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.fna.gz" \
    "Variant1" "https://example.com/DENV_2_comparison1.fa"
process_comparison_genome "DENV-2" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/871/845/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.fna.gz" \
    "Variant2" "https://example.com/DENV_2_comparison2.fa"
process_comparison_genome "DENV-2" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/871/845/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.fna.gz" \
    "Variant3" "https://example.com/DENV_2_comparison3.fa"

# Process DENV-3 Genomes
process_reference_genome "DENV-3" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788295.1_ASM478829v1/GCF_004788295.1_ASM478829v1_genomic.fna.gz"
process_reference_annotation "DENV-3" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788295.1_ASM478829v1/GCF_004788295.1_ASM478829v1_genomic.gff.gz"

process_comparison_genome "DENV-3" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788295.1_ASM478829v1/GCF_004788295.1_ASM478829v1_genomic.fna.gz" \
    "Variant1" "https://example.com/DENV_3_comparison1.fa"
process_comparison_genome "DENV-3" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788295.1_ASM478829v1/GCF_004788295.1_ASM478829v1_genomic.fna.gz" \
    "Variant2" "https://example.com/DENV_3_comparison2.fa"
process_comparison_genome "DENV-3" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788295.1_ASM478829v1/GCF_004788295.1_ASM478829v1_genomic.fna.gz" \
    "Variant3" "https://example.com/DENV_3_comparison3.fa"

# Process DENV-4 Genomes
process_reference_genome "DENV-4" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788296.1_ASM478830v1/GCF_004788296.1_ASM478830v1_genomic.fna.gz"
process_reference_annotation "DENV-4" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788296.1_ASM478830v1/GCF_004788296.1_ASM478830v1_genomic.gff.gz"

process_comparison_genome "DENV-4" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788296.1_ASM478830v1/GCF_004788296.1_ASM478830v1_genomic.fna.gz" \
    "Variant1" "https://example.com/DENV_4_comparison1.fa"
process_comparison_genome "DENV-4" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788296.1_ASM478830v1/GCF_004788296.1_ASM478830v1_genomic.fna.gz" \
    "Variant2" "https://example.com/DENV_4_comparison2.fa"
process_comparison_genome "DENV-4" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788296.1_ASM478830v1/GCF_004788296.1_ASM478830v1_genomic.fna.gz" \
    "Variant3" "https://example.com/DENV_4_comparison3.fa"

echo "All reference genomes, annotations, and comparison genomes have been processed and added to JBrowse."
