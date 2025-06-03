Analyze allele and genotype frequencies of tetraploid SNP data
with optional binning of observed frequencies for smoother curves.

This script reads samples from a VCF file and computes allele frequency 
vs. genotype frequency distributions for selected samples. It outputs
a plot of genotype frequency vs. allele frequency in SVG file format, showing:

- rasterized observed datapoints,
- raw and smoothed mean data curves, 
- theoretical curves for autotetraploid (tetrasomic) and allotetraploid (disomic) models.

usage:
- tetraploid_analyses.py [-h] -s SAMPLES [SAMPLES ...] [-o OUTPUT]
                              [--binsize BINSIZE] [--dpi DPI] [--smooth]
                              [--smooth-factor SMOOTH_FACTOR]
                              vcf_file
