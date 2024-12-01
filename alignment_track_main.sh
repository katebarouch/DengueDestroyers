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

# Function to process a reference genome
process_reference_genome() {
    virus_name="$1"
    genome_url="$2"
    assembly_name="$3" 

    echo "Processing $virus_name reference genome..."
    wget "$genome_url" -O "$WORKDIR/${virus_name}_genome.fasta"
    check_error

    echo "Indexing $virus_name reference genome..."
    samtools faidx "$WORKDIR/${virus_name}_genome.fasta"
    check_error

    echo "Adding $virus_name reference genome to JBrowse with assembly name: $assembly_name..."
    jbrowse add-assembly "$WORKDIR/${virus_name}_genome.fasta" --out "$APACHE_ROOT/jbrowse2" --load copy --name "$assembly_name"
    check_error
}

# Function to process a reference annotation
# needs to be fixed if we no longer recieve zipped files... also don't have links yet
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

# Function to process comparison genome for alignment track
process_comparison_genome() {
    reference_virus_name="$1"
    comparison_virus_name="$2"
    comparison_genome_url="$3"
    assembly_name="$4"

    echo "Building Bowtie2 index for $reference_virus_name genome..."
    bowtie2-build "$WORKDIR/${reference_virus_name}_genome.fasta" "$WORKDIR/${reference_virus_name}_genome_index"
    check_error

    echo "Downloading $comparison_virus_name..."
    wget -q "$comparison_genome_url" -O "$WORKDIR/${comparison_virus_name}.fasta"
    check_error

    echo "Aligning $comparison_virus_name to $reference_virus_name reference genome..."
    bowtie2 -x "$WORKDIR/${reference_virus_name}_genome_index" -f "$WORKDIR/${comparison_virus_name}.fasta" -S "$WORKDIR/${comparison_virus_name}.sam"
    check_error

    echo "Converting $comparison_virus_name alignment to BAM format..."
    samtools view -bS "$WORKDIR/${comparison_virus_name}.sam" > "$WORKDIR/${comparison_virus_name}.bam"
    check_error

    echo "Sorting $comparison_virus_name BAM file..."
    samtools sort "$WORKDIR/${comparison_virus_name}.bam" -o "$WORKDIR/${comparison_virus_name}.sorted.bam"
    check_error

    echo "Indexing $comparison_virus_name sorted BAM file..."
    samtools index "$WORKDIR/${comparison_virus_name}.sorted.bam"
    check_error

    echo "Adding $comparison_virus_name alignment track to JBrowse..."
    jbrowse add-track "$WORKDIR/${comparison_virus_name}.sorted.bam" --out "$APACHE_ROOT/jbrowse2" --load copy --name "$comparison_virus_name" --assemblyNames "$assembly_name"
    check_error
}

### Main Script Execution ###

# Process DENV-1 Genomes
process_reference_genome "DENV-1_New_Caledonia-2017-AVS-NC-094" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW315194.1&report=fasta&format=text" "DENV-1"
# process_reference_annotation "DENV-1" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz"

process_comparison_genome "DENV-1_New_Caledonia-2017-AVS-NC-094" \
    "DENV-1_Xishuangbanna_Dai_China" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW386867.1&report=fasta&format=text" "DENV-1"
process_comparison_genome "DENV-1_New_Caledonia-2017-AVS-NC-094" \
    "DENV-1_Guangzhou_China" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PQ357572.1&report=fasta&format=text" "DENV-1"

# Process DENV-2 Genomes
process_reference_genome "DENV-2_Cosmopolitan" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW512496.1&report=fasta&format=text" "DENV-2"
# process_reference_annotation "DENV-2, Cosmopolitan" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/871/845/GCF_000871845.1_ViralProj20183/GCF_000871845.1_ViralProj20183_genomic.gff.gz"

process_comparison_genome "DENV-2_Cosmopolitan" \
    "DENV-2_Malaysia" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PQ465587.1&report=fasta&format=text" "DENV-2"
process_comparison_genome "DENV-2_Cosmopolitan" \
    "DENV-2_Brazil" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=PP546320.1&report=fasta&format=text" "DENV-2"

# Process DENV-3 Genomes
process_reference_genome "DENV-3_Florida" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OQ821613.1&report=fasta&format=text" "DENV-3"
# process_reference_annotation "DENV-3" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788295.1_ASM478829v1/GCF_004788295.1_ASM478829v1_genomic.gff.gz"

process_comparison_genome "DENV-3_Florida" \
    "DENV-3_V91A_mutation" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OQ821525.1&report=fasta&format=text" "DENV-3"
process_comparison_genome "DENV-3_Florida" \
    "DENV-3_Senegal" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=MW288026.1&report=fasta&format=text" "DENV-3"

# Process DENV-4 Genomes
process_reference_genome "DENV-4_Indonesia" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OL314742.1&report=fasta&format=text" "DENV-4"
# process_reference_annotation "DENV-4, Indonesia" "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/788/295/GCF_004788296.1_ASM478830v1/GCF_004788296.1_ASM478830v1_genomic.gff.gz"

process_comparison_genome "DENV-4_Indonesia" \
    "DENV-4_Paraguay" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OP811943.1&report=fasta&format=text" "DENV-4"
process_comparison_genome "DENV-4_Indonesia" \
    "DENV-4_Cambodia_2010" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=KF543272.1&report=fasta&format=text" "DENV-4"
process_comparison_genome "DENV-4_Indonesia" \
    "DENV-4_French_Polynesia_2009" "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=JN832541.1&report=fasta&format=text" \ "DENV-4"
    

echo "All reference genomes, annotations, and comparison genomes have been processed and added to JBrowse."
