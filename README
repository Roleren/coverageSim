# coverageSim
### Simulating Ribosome profiling data for tool validation
![](inst/images/coverageSim_overview.png)


## How to install
```r
library(devtools)
install_github("Roleren/coverageSim")
```

## Simple example

Lets make a genome with 6 genes, each being coding (they have a CDS) and
having 1 translated uORF:

```r
library(coverageSim)
library(ORFik)
## Simple example
# 6 genes on 6 chromosomes (1 active uORF each)
simGenome6 <- simGenome(n = 6, max_uorfs = 1)
# Simulate Ribo-seq only
gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
libtypes = "RFP", print_statistics = FALSE)
region_count_table <- simCountTablesRegions(gene_count_table,
     regionsToSample = c("leader", "cds", "trailer"))
simNGScoverage(simGenome6, region_count_table)
```
