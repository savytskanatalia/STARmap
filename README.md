# STARmap
Simple nextflow pipeline to map reads and, if needed, generate genome index


Example command:

```
nextflow run scripts/star_map.nf -c scripts/star_map.config \\
         --reads "/my/full/dir/*.fastq.gz" \\
         -with-docker savytskanatalia/quantefication2 \\
         --star_index "/my/full/dir/star-index/" \\
         --oufiltmm 1 --winamm 1 --singleEnd
```
If supplied with Star index, runs mapping directly. If supplied with genome gtf+fasta starts with generating genome index first.
