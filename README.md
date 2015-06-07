# BAMreport
Report on a BAM file. 

BAMreport creates a directory (which can be user specified) and puts all the 
output files in it. The output files include the following:

* stats.txt : Alignment statistics
* nucleotides.pdf : Nucleotide composition vs. position on sequence
* quality.pdf : Base quality vs position on sequence
* coverage.pdf : Depth of coverage distribution
* gccoverage.png : Variation of GC and coverage 
* insertlengths.pdf : Insert length distribution
* readlengths.pdf : Variation of read lengths

## REQUIREMENTS
These tools should work on any standard 64 bit Linux environment with
* GCC
* Python (version >= 2.7.6)

## TEST DATA
The test dataset includes a BAM file and a fasta reference file. Run the tool 
on the dataset by running the following:

```
    ./../src/BAMreport -w 100 reference.fa alignments.bam
```

This will create an output directory "reports" with all the files in it. These
files should be the same as the ones included in the "expected" directory.
