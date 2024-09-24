# Description
This Nextflow pipeline extracts metadata from GEO series and run info from SRA and join them together.

# example
```
nextflow run main.nf --sra_id 'SRP161855' --geo_series 'GSE120011' -with-singularity -resume
```
