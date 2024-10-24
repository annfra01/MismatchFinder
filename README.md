# MismatchFinder
A tool to find variants (mismatched base pairs, insertions, deletions) in genomic sequences.

## Getting started 
These instructions will help you copy the project to run it on your local machine.

### Prerequisites

The source code can be cloned to your local directory using:
```
git clone https://github.com/annfra01/MismatchFinder.git
```
It can also be downloaded and extracted from the github page: https://github.com/annfra01/MismatchFinder

The file requirements.txt can be used to install all necessary dependencies. It is recommended to setup a virtual environment for this project.

Once you have setup your virtual environment, run the following code to install the dependencies:
```
pip install -r requirements.txt
```
or if your pip points to an existing Python2 environment:
```
pip3 install -r requirements.txt
```


## How to use the tool
You can start the tool by using this command:
```
python3 MismatchFinder.py -f example.bam -i example.bam.bai
```

### Input parameters
+ -f : A .bam file.
+ -i : A .bai file. It should be generated beforehand using the following command:

```
samtools index example.bam example.bam.bai
```
### Output
The tool will output:
1. a .bed file mismatches.bed with 6 columns. The columns are: chromosome, start position of mismatch, end position of mismatch, type of mismatch (regular mismatch, insertion, deletion), a dot as a placeholder, the direction of the read (- is reverse, + is forward). The .bed output file for example.bam looks like this:
   
   <img src="/example_bed.png" height="150">

3. statistics using the standard output. The statistics for the example.bam looks like this:

   <img src="/example_statistics.png" height="150">

   
