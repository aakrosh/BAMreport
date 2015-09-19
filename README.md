# BAMreport
Report on a BAM file. 

BAMreport creates a PDF report (which can be user specified) and puts all the 
QC graphs and plots in it. The output file includes the following:

* Alignment statistics from the BAM file. 
* Nucleotide composition on both reads of fragment as the position on sequence
  is varied.
* Base quality variation with position on sequence.
* Depth of coverage distribution.
* Relative depth of coverage distribution based on the specified window size.
* Variation of GC and coverage. 
* Insert length distribution.
* Variation of read lengths.
* Frequency of clipping on the 5' end which is indicative of adapters or other
  contaminations.
* Mismatches on the reads against the reference and their variation with
  position.
* Indels on the reads against the reference and their variation with position on
  the reads.

## REQUIREMENTS
These tools should work on any standard 64 bit Linux environment with
* GCC
* Python (version >= 2.7.6)

## TEST DATA
The test dataset includes a BAM file and a fasta reference file. Run the tool 
on the dataset by running the following:

```
    ./../src/BAMreport -w 100  reference.fa alignments.bam
```

This should create a file "report.pdf" which should have the output similar to 
expected.pdf in the test_data directory.
