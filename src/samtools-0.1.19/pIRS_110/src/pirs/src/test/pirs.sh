./../pirs diploid -i ref_seq.fa -s 0.001 -a 2 -d 0.0001 -v 0.000001 -c 1 -o ref22 >SimDiploid.out 2>SimDiploid.err
./../pirs simulate -i ref_seq.fa  -M 1 -m 800 -l 100 -x 3 -v 40 >SimReads.out 2>SimReads.err
./../pirs simulate -E 1 -o EAMSS2 -e 0.01 -i ref_seq.fa -I ref22.snp.indel.inversion.fa.gz -M 1 -m 800 -l 90 -x 3 -v 40 >SimReadsEAMSS2.out 2>SimReadsEAMSS2.err

