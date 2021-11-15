
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


## Setting Up R Environment


``` {R environment, message = FALSE}
library(dada2, quietly = TRUE)
packageVersion("dada2")
```


``` {R filenames}
path <- "reads" #set this to the path where your fastq files live
list.files(path)
forward_reads <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
reverse_reads <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
samples <- sapply(strsplit(basename(forward_reads), "_"), `[`, 1)
```


## Quality Plot Inspection


```{r plotQuality}
#currently only running the first 4 of the list to save time
plotQualityProfile(forward_reads[1:4])
plotQualityProfile(reverse_reads[1:4])
```


## Filter and Trimming


```{r filterNames}
filtered_forward_reads <- paste0(samples, "_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fq.gz")
```

```{r filterAndTrim}
filtered_out <- filterAndTrim(forward_reads, 
                              filtered_forward_reads,
                              reverse_reads, 
                              filtered_reverse_reads, 
                              maxEE=c(2,2),
                              minLen=175, 
                              truncLen=c(250,200))
```

```{r ViewFiltered}
filtered_out
```

```{r plotQualityFiltered}
plotQualityProfile(filtered_forward_reads[1:4])
plotQualityProfile(filtered_reverse_reads[1:4])
```


## Generate Error Model


```{r errorModel}
err_forward_reads <- learnErrors(filtered_forward_reads)
err_reverse_reads <- learnErrors(filtered_reverse_reads)
#set multithread = TRUE if running on your own system:
#err_forward_reads <- learnErrors(filtered_forward_reads, multithread = TRUE)
#err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread = TRUE)
```

```{r plotErrors}
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
```

## Dereplication

```{r dereplication}
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples 
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
```

# Inferring ASVs

```{r inferASV}
#dada_forward <- dada(derep_forward, err=err_forward_reads, #pool="pseudo")
#dada_reverse <- dada(derep_reverse, err=err_reverse_reads, #pool="pseudo")
#set multithread = TRUE if running on your own system:
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread = TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread = TRUE)
```


## Merging paired reads

```{r merge}
merged_amplicons <- mergePairs(dada_forward, 
                              derep_forward, 
                              dada_reverse,
                              derep_reverse, 
                              minOverlap=20)
```


## Count Table and Summary

```{r seqtab}
seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
```

```{r removeChimeras}
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

```{r readCounts}
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
               filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
               dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
               nonchim=rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)
```

