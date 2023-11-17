# MismatchFinder
A tool to find variants (regular mismatches, Insertions, Deletions) in genomic sequences.
## How to use the tool
You can start the tool by using this command:
```
MismatchFinder.py -f file.bam -i file.bam.bai
```

## Input parameters
+ -f : A .bam file.
+ -i : A .bai file. It should be generated beforehand using the following command:

```
samtools index file.bam file.bam.bai
```
## Output
The tool will output:
1. a .bed file which looks like this:
2. a statistics using the standard output which looks like this:
   
