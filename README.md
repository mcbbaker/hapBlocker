# hapBlocker
**_(Work Related)_**

Algorithm to output haplotype blocks across a population in the manner of a mosaic-style plot. 

## Overview
Given a VCF with samples from a population, it is desirable to find the genomic regions (blocks) where there is a high percentage of the same SNPs between samples. This algorithm takes in a VCF from a population and outputs these blocks in the manner of a mosaic-style plot. 

This 'mosaic-style' means that the first sample will have one block with it's own haplotype for each chromosome. For the remaining samples, the algorithm compares the SNPs from the current sample to the above samples and if there are enough identical SNPs (and other criteria passes), the current sample will declare a block with the above sample's haplotype. Precedence for comparison is given to the top sample, working downward. If no above sample passes the criteria, a block is declared with the sample's haplotype. 

## Usage & Parameters
python3 hapBlocker.py vcf.gz index.gz.[csi,tbi] samples min_snps percent-threshold comparators

**_Parameters (all required):_**
- **_vcf.gz:_** gzipped VCF file. Note that chromosomes must be named such that 'Chr' is included
- **_index.gz:_** gzipped VCF index file in csi or tbi format
- **_samples:_** path to a file that lists the samples, each on their own line. The names of the samples must match the names in the VCF file. If wanting to visualize, list the samples in the order they should appear in the plot
- **_min_snps:_** the number of SNPs that must be compared before a block can be declared
- **_percent-threshold:_** the percent of identical SNPs across 'min_snps' needed to declare a block
- **_comparators:_** the first X rows in the 'samples' file that should be compared to. For example, if you have many samples and only want to see blocks with parental haplotypes, choose 2 (assuming parents are the first 2 samples in the 'samples' file). If you want samples to be compared to every sample above it, choose the number of samples-1

## Output
The blocks are outputted to stdout in the following format:
 > Chr  Start  End  Block#  Sample  Haplotype

The output is tab-delimited, where 'Block#' refers to the nth outputted block for the sample 
 
