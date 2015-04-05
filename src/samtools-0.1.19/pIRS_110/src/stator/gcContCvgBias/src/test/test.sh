./../gc_coverage_bias -r gcdeptestref.fa.gz -o test -w 100,200 --gcdump --depwindump gcdeptest.depth.gz
perl ./../gc_coverage_bias_plot.pl test_100.dat
perl ./../gc_coverage_bias_plot.pl test_200.dat
ls -l test*
