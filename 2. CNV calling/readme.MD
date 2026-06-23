Here are the CNV calling scripts located:

1. CNVkit script takes batch or single BAM as input.

2. It looks for reference genome fasta file and blacklist coordinates file, if it doesn't find them, it downloads them.

3. It optionally looks for pooled reference files.

4. It performs CNV calling and segmentation with user-defined method.

5. It isolates background low coverage regions and sends them to be processed by QDNAseq.

6. It combines the results together.
