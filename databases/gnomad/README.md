#gnomAD
gnomAD v3 released! 71,702 genomes aligned on GRCh38.
## download
```
make all -j 24
```

## Filter variants
there are variants that do not pass all the filter, we have to keept only those passing filters.
## Exome variants

## Basic Stats of passing filter variants

### Number of SNPs
526,001,545 (526 million)
```
grep "number of SNPs:" *.stats  | awk '{i+=$NF;}END{print i}'
```
### Number of Indels
69,168,024 (69 million)
```
grep "number of indels:" *.stats  | awk '{i+=$NF;}END{print i}'
```
### Total variants
595,169,569 millions variantss.
