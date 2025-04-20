# Demos - veupathdb


## convertSnpCacheToVCF.pl [cachedSnpFile] [referencStrain]

Writes VCF to stdout!

Example:
```
convertSnpCacheToVCF.pl examples/isolateSNPs.cache 3D7 > output.vcf

Reformat a tab-delimited file of SNPs (from arrays) into a VCF file.
A vcf file could be made for each strain OR one file for the entire dataset.  For our purpose we'd want one file for the entire dataset with columns for each strain.  

[VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

TODO:
* Support gzipped input file
* Use GetOptions (GetOpt::Long) instead of ordered arguments
* Support diploid alleles
* Add more INFO data fields
* Create a generic writer package
