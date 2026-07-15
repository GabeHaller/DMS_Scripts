#!/bin/bash

# Set the path to bbmerge
BBMERGE=/neuro/haller/hallerg/Meade/Meade/bbmap/bbmerge.sh

# Loop over all R1 FASTQ files
for r1 in *_R1_001.fastq.gz; do
    # Get base name (remove _R1_001.fastq.gz)
    base=${r1%%_R1_001.fastq.gz}
    
    # Define matching R2 file
    r2="${base}_R2_001.fastq.gz"
    
    # Output file name
    out="${base}_Merged.fastq.gz"
    
    # Run bbmerge
    echo "Merging $r1 and $r2 → $out"
    "$BBMERGE" in1="$r1" in2="$r2" out="$out"
done
