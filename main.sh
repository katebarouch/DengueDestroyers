#!/bin/bash

### Dengue ###

#!/bin/bash

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

### Download and Process the Reference Genome ###
echo "Downloading Dengue genome (reference genome)..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.fna.gz" -O $WORKDIR/genome.fna.gz
check_error

echo "Decompressing reference genome file..."
gunzip -f $WORKDIR/genome.fna.gz
check_error

echo "Renaming reference genome file..."
mv $WORKDIR/genome.fna $WORKDIR/viral_genome.fa
check_error

echo "Indexing reference genome file with samtools..."
samtools faidx $WORKDIR/viral_genome.fa
check_error

echo "Adding single-stranded RNA reference genome assembly to JBrowse..."
jbrowse add-assembly $WORKDIR/viral_genome.fa --out $APACHE_ROOT/jbrowse2 --load copy
check_error

### Download and Process Annotations ###
echo "Downloading genome annotations..."
wget -q "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/862/125/GCF_000862125.1_ViralProj15306/GCF_000862125.1_ViralProj15306_genomic.gff.gz" -O $WORKDIR/annotations.gff.gz
check_error

echo "Decompressing annotations file..."
gunzip -f $WORKDIR/annotations.gff.gz
check_error

echo "Sorting GFF3 annotations..."
jbrowse sort-gff $WORKDIR/annotations.gff > $WORKDIR/genes.gff
check_error

echo "Compressing sorted GFF3 file..."
bgzip -f $WORKDIR/genes.gff
check_error

echo "Indexing compressed GFF3 file with tabix..."
tabix $WORKDIR/genes.gff.gz
check_error

echo "Adding annotations track to JBrowse..."
jbrowse add-track $WORKDIR/genes.gff.gz --out $APACHE_ROOT/jbrowse2 --load copy
check_error

echo "Reference dengue genome and annotations successfully added to JBrowse."

### Align Comparison Genome to Reference Genome ###

echo "Preparing reference genome..."
bowtie2-build $WORKDIR/viral_genome.fa $WORKDIR/viral_genome
check_error

echo "Downloading and decompressing comparison genome..."
wget "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=OR486055.1&report=fasta&format=text" -O dengue_virus1.fa
check_error 

echo "Aligning comparison genome to reference genome using Bowtie2..."
bowtie2 -x $WORKDIR/viral_genome -f $WORKDIR/dengue_virus1.fa -S $WORKDIR/comparison_genome.sam
check_error 

echo "Converting SAM file to BAM file..."
samtools view -bS $WORKDIR/comparison_genome.sam > $WORKDIR/comparison_genome.bam
check_error 

echo "Sorting BAM file..."
samtools sort $WORKDIR/comparison_genome.bam -o $WORKDIR/comparison_genome.sorted.bam
check_error 

echo "Indexing BAM file..."
samtools index $WORKDIR/comparison_genome.sorted.bam
check_error 

echo "Adding alignment track to JBrowse..."
jbrowse add-track $WORKDIR/comparison_genome.sorted.bam --out $APACHE_ROOT/jbrowse2 --load copy
check_error 

echo "Comparison genome alignment successfully added to JBrowse."
