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

# Define URLs and JBrowse names for comparison genomes (all 3 entries for each serotype are stored in single string for simplicity)
declare -A comparison_genomes
comparison_genomes["DENV-1"]="https://example.com/DENV_1_comparison1.fa,SampleA https://example.com/DENV_1_comparison2.fa,SampleB https://example.com/DENV_1_comparison3.fa,SampleC"
comparison_genomes["DENV-2"]="https://example.com/DENV_2_comparison1.fa,IsolateX https://example.com/DENV_2_comparison2.fa,IsolateY https://example.com/DENV_2_comparison3.fa,IsolateZ"
comparison_genomes["DENV-3"]="https://example.com/DENV_3_comparison1.fa,Variant1 https://example.com/DENV_3_comparison2.fa,Variant2 https://example.com/DENV_3_comparison3.fa,Variant3"
comparison_genomes["DENV-4"]="https://example.com/DENV_4_comparison1.fa,StrainAlpha https://example.com/DENV_4_comparison2.fa,StrainBeta https://example.com/DENV_4_comparison3.fa,StrainGamma"

# Loop over each reference genome and its corresponding annotations
for virus in DENV-1 DENV-2 DENV-3 DENV-4; do
  # Process reference genome
  echo "Processing $virus reference genome..."
  wget -q "${genome_urls[$virus]}" -O "$WORKDIR/${virus}_genome.fna.gz"
  check_error

  echo "Unzipping $virus reference genome..."
  gunzip -f "$WORKDIR/${virus}_genome.fna.gz"
  check_error

  echo "Indexing $virus reference genome..."
  samtools faidx "$WORKDIR/${virus}_genome.fna"
  check_error
  
  echo "Adding $virus reference genome to JBrowse..."
  jbrowse add-assembly "$WORKDIR/${virus}_genome.fna" --out "$APACHE_ROOT/jbrowse2" --load copy
  check_error

  # Process reference annotations
  echo "Processing $virus annotations..."
  wget -q "${annotation_urls[$virus]}" -O "$WORKDIR/${virus}_annotations.gff.gz"
  check_error

  echo "Unzipping $virus annotations..."
  gunzip -f "$WORKDIR/${virus}_annotations.gff.gz"
  check_error

  echo "Sorting $virus annotations..."
  jbrowse sort-gff "$WORKDIR/${virus}_annotations.gff" > "$WORKDIR/${virus}_genes.gff"
  check_error

  echo "Compressing $virus annotations..."
  bgzip -f "$WORKDIR/${virus}_genes.gff"
  check_error

  echo "Indexing $virus annotations..."
  tabix "$WORKDIR/${virus}_genes.gff.gz"
  check_error

  echo "Adding $virus annotations to JBrowse..."
  jbrowse add-track "$WORKDIR/${virus}_genes.gff.gz" --out "$APACHE_ROOT/jbrowse2" --load copy
  check_error

  echo "$virus genome and annotations successfully added to JBrowse."

  # Align comparison genomes
  echo "Building Bowtie2 index for $virus reference genome..."
  bowtie2-build "$WORKDIR/${virus}_genome.fna" "$WORKDIR/${virus}_genome_index"
  check_error

  # Process each comparison genome
  echo "Processing comparison genomes for $virus..."
  for comparison_entry in $(echo "${comparison_genomes[$virus]}"); do
    comparison_url=$(echo "$comparison_entry" | cut -d',' -f1)
    jbrowse_name=$(echo "$comparison_entry" | cut -d',' -f2)

    echo "Downloading $jbrowse_name for $virus..."
    wget -q "$comparison_url" -O "$WORKDIR/${jbrowse_name}.fa"
    check_error

    echo "Aligning $jbrowse_name to $virus reference genome..."
    bowtie2 -x "$WORKDIR/${virus}_genome_index" -f "$WORKDIR/${jbrowse_name}.fa" -S "$WORKDIR/${jbrowse_name}.sam"
    check_error

    echo "Converting $jbrowse_name alignment to BAM format..."
    samtools view -bS "$WORKDIR/${jbrowse_name}.sam" > "$WORKDIR/${jbrowse_name}.bam"
    check_error

    echo "Sorting $jbrowse_name BAM file..."
    samtools sort "$WORKDIR/${jbrowse_name}.bam" -o "$WORKDIR/${jbrowse_name}.sorted.bam"
    check_error

    echo "Indexing $jbrowse_name sorted BAM file..."
    samtools index "$WORKDIR/${jbrowse_name}.sorted.bam"
    check_error

    echo "Adding $jbrowse_name alignment track to JBrowse..."
    jbrowse add-track "$WORKDIR/${jbrowse_name}.sorted.bam" --out "$APACHE_ROOT/jbrowse2" --load copy --name "$jbrowse_name"
    check_error
  done
  echo "All comparison genomes for $virus successfully aligned and added to JBrowse."
done
