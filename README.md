# DMS_Scripts

These are the scripts used in the SGCA DMS paper. 
1. Merge FASTQ files with BBMerge using Merge_BBmerge_Fastq. The location of the BBMerge program will need to be updated in the .sh file. This requires that the R1 and R2 read files are overlapping sufficiently to merge the reads into one large file. I.e. the insert is less than 2x the length of the reads, for example 250bp insert for 2x150bp paired-end reads.
2. Align using BWA to a cDNA Fasta file that starts with the ATG start site for the gene bring mutated.
3. Using the .sam files produced by BWA, run Count_AAchanges_SamFasta_CodonCounts2.pl using the cDNA Fasta and Merged FASTQ files as input.
4. You can then merge multiple count files using the Concatenate_CountLists.pl program. 
