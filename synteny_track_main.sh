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
process_genome() {
    virus_name="$1"
    genome_url="$2"

    echo "Processing $virus_name reference genome..."
    wget "$genome_url" -O "$WORKDIR/${virus_name}_genome.fasta"
    check_error

    echo "Indexing $virus_name reference genome..."
    samtools faidx "$WORKDIR/${virus_name}_genome.fasta"
    check_error

    echo "Adding $virus_name reference genome to JBrowse..."
    jbrowse add-assembly "$WORKDIR/${virus_name}_genome.fasta" --out "$APACHE_ROOT/jbrowse2" --load copy --name "$virus_name"
    check_error
}

# Function to process a reference annotation
#needs to be fixed
process_annotation() {
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

# Function to process comparison genome for synteny track
create_synteny_track() {
    reference_virus_name="$1"
    comparison_virus_name="$2"

    echo "Creating synteny track between $reference_virus_name and $comparison_virus_name..."
    minimap2 "$WORKDIR/${reference_virus_name}_genome.fasta" "$WORKDIR/${comparison_virus_name}_genome.fasta" > "$WORKDIR/${reference_virus_name}_vs_${comparison_virus_name}.paf"
    check_error

    echo "Adding synteny track to JBrowse..."
    jbrowse add-track "$WORKDIR/${reference_virus_name}_vs_${comparison_virus_name}.paf" \
        --assemblyNames "$reference_virus_name,$comparison_virus_name" \
        --load copy \
        --out "$APACHE_ROOT/jbrowse2"
    check_error
}

### Main Script Execution ###

# Process DENV-1 Genomes
process_genome "DENV-1_New_Caledonia-2017-AVS-NC-094" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW315194.1&report=fasta&format=text" 
process_genome "DENV-1_Xishuangbanna_Dai_China" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW386867.1&report=fasta&format=text" 
process_genome "DENV-1_Guangzhou_China" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PQ357572.1&report=fasta&format=text" 


create_synteny_track "DENV-1_New_Caledonia-2017-AVS-NC-094" \
    "DENV-1_Xishuangbanna_Dai_China" 
create_synteny_track "DENV-1_Xishuangbanna_Dai_China" \
    "DENV-1_Guangzhou_China" 
create_synteny_track "DENV-1_Guangzhou_China" \
    "DENV-1_New_Caledonia-2017-AVS-NC-094" 

# Process DENV-2 Genomes
process_genome "DENV-2_Cosmopolitan" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW512496.1&report=fasta&format=text" 
process_genome "DENV-2_Malaysia" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PQ465587.1&report=fasta&format=text" 
process_genome "DENV-2_Brazil" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PP546320.1&report=fasta&format=text" 

create_synteny_track "DENV-2_Cosmopolitan" \
    "DENV-2_Malaysia" 
create_synteny_track "DENV-2_Cosmopolitan" \
    "DENV-2_Brazil" 
create_synteny_track "DENV-2_Malaysia" \
    "DENV-2_Brazil" 

# Process DENV-3 Genomes
process_genome "DENV-3_Florida" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OQ821613.1&report=fasta&format=text" 
process_genome "DENV-3_V91A_mutation" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OQ821525.1&report=fasta&format=text" 
process_genome "DENV-3_Senegal" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW288026.1&report=fasta&format=text" 

create_synteny_track "DENV-3_Florida" \
    "DENV-3_V91A_mutation" 
create_synteny_track "DENV-3_Florida" \
    "DENV-3_Senegal" 
create_synteny_track "DENV-3_V91A_mutation" \
    "DENV-3_Senegal" 

# Process DENV-4 Genomes
process_genome "DENV-4_Indonesia" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OL314742.1&report=fasta&format=text" 
process_genome "DENV-4_Paraguay" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OP811943.1&report=fasta&format=text" 
process_genome "DENV-4_Cambodia_2010" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=KF543272.1&report=fasta&format=text" 
process_genome "DENV-4_French_Polynesia_2009" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=JN832541.1&report=fasta&format=text" 

create_synteny_track "DENV-4_Indonesia" \
    "DENV-4_Paraguay" 
create_synteny_track "DENV-4_Indonesia" \
    "DENV-4_Cambodia_2010" 
create_synteny_track "DENV-4_Indonesia" \
    "DENV-4_French_Polynesia_2009" 
create_synteny_track "DENV-4_Paraguay" \
    "DENV-4_Cambodia_2010" 
create_synteny_track "DENV-4_Cambodia_2010" \
    "DENV-4_French_Polynesia_2009" 

echo "All reference genomes, annotations, and comparison genomes have been processed and added to JBrowse."
