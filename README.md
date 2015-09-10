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

Here I show the expected output of the test dataset and explain each of them briefly
* stats.txt 
```
Reads Summary
-------------
QC-passed reads : 2630
QC-failed reads : 0
Number of single-fragments           : 0
Number of reads paired in sequencing : 2630

Duplicate reads         : 0 (0.00%)
Number of reads aligned : 2630 (100.00%)
Number of reads mapped  : 2630 (100.00%)

Number of read1's : 1315
Number of read2's : 1315
Number of reads properly paired              : 2630
Number of reads with itself and mate aligned : 2630

Bases Summary
-------------
QC-passed bases : 199880
QC-failed bases : 0
Number of bases in single-fragments  : 0
Number of bases paired in sequencing : 199880

Duplicate bases         : 0 (0.00%)
Number of bases aligned : 199880 (100.00%)
Number of bases mapped  : 199880 (100.00%)
```
This file summarizes the alignment statistics. "Mapped" refers to the reads that align uniquely i.e. with a mapping quality > 0. The aligned and mapped count includes the putative PCR duplicates.
* readlengths.pdf : This file shows the distribution of readlengths in your BAM file. The left-panel shows read1's and the right panel shows the length distribution of read2's.
* nucleotides.pdf : This file shows the nucleotide composition of the sequences vs. the position on the sequence. Again, the left-panel shows read1's and the right panel shows the read2's. Ideally you do not want to see any bias in the nucleotide composition with position on the sequences.
* quality.pdf : This file shows the distribution of the quality values with respect to the position on the sequences. We show the first, second, and the third quantile values as well as the minimum and maximum observed quality values at every location.
* insertlengths.pdf : This file shows the distribution of insert length (outer distance) between the paired-end sequences.
* coverage.pdf : This shows the depth-of-coverage distribution of the sequences.
* gccoverage.pdf : This shows a scatterplot of average coverage vs. the average GC content in non-overlapping windows on the reference. 
* fiveprimeclips.pdf : This shows a histogram of the clipping location on the 5'
 end of the sequences. Clipping on the 5' end is frequently a signature of
adapter sequences or some other contamination. Clipping on the 3' end on the
other hand can happen due to trimming of low quality regions.
