# BAMreport
Report on a BAM file. 

BAMreport creates a directory (which can be user specified) and puts all the 
output files in it. The output files include the following:

1) stats.txt : Alignment statistics
2) nucleotides.pdf : Nucleotide composition vs. position on sequence
3) quality.pdf : Base quality vs position on sequence
4) coverage.pdf : Depth of coverage distribution
5) gccoverage.png : Variation of GC and coverage 
6) insertlengths.pdf : Insert length distribution
7) readlengths.pdf : Variation of read lengths

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
