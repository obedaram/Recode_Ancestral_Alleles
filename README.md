# This script recodes your VCF file by ancestral allele.
The Ancestral allele files provided were queried from the 1000 genomes AA field in the VCF files. These sites are all biallelic SNPs and only contains sites that were determined that it can be successfully flipped, i.e. the AA field matches either the REF or the ALT.
This script will automatically filter and only keep those sites intersecting the Ancestral Allele text file and print out a summary of what was done as well as a file of how each allele was treated. The script will accept gzipped or normal txt files for the Ancestral Alleles, and will accept uncompressed VCF or bgzip compressed VCF
Usage:
```
python recodeAA.py <ancestral_allele.txt> <vcf_file.gz>
```
