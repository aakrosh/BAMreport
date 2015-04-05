all: configure simulator stator symlinks

PERL_LIST = stator/baseCallingMatrix/baseCalling_Matrix_analyzer \
	stator/baseCallingMatrix/baseCalling_Matrix_calculator \
	stator/baseCallingMatrix/baseCalling_Matrix_merger \
	stator/alignment_stator \
	stator/indelstat_sam_bam \
	stator/gcContCvgBias/src/gc_coverage_bias_plot

configure:
	-mkdir -p ${PREFIX} 2>/dev/null
	cd ./src/pirs/ && ./configure --prefix=${PREFIX}
	cd ./src/stator/gcContCvgBias/ && ./configure --prefix=${PREFIX}

install: symlinks
	-mkdir -p ${PREFIX}
	find . -type l|while read a;do cp --copy-contents -LR  "$$a" ${PREFIX}/;done

simulator:
	cd ./src/pirs/ && ${MAKE}
  
stator:
	cd ./src/stator/gcContCvgBias/ && ${MAKE}

symlinks:
	-@ln -s ./src/pirs/src/pirs 2> /dev/null
	-@ln -s ./src/pirs/src/Profiles 2> /dev/null
	-@ln -s ./src/stator/gcContCvgBias/src/gc_coverage_bias 2> /dev/null
	-@ln -s ./src/stator/alignment_stator.pl alignment_stator 2> /dev/null
	-@for P in ${PERL_LIST}; do \
	    ln -s ./src/$${P}.pl `basename $${P}` 2> /dev/null; \
	done

test: 
	cd ./src/pirs/src/test && sh pirs.sh
	cd ./src/stator/gcContCvgBias/src/test && sh test.sh

clean:
	cd ./src/pirs && ${MAKE} clean
	cd ./src/pirs/src && ${MAKE} clean
	cd ./src/stator/gcContCvgBias && ${MAKE} clean
	cd ./src/stator/gcContCvgBias/src && ${MAKE} clean
	-rm pirs gc_coverage_bias Profiles

allclean: clean
	cd ./src/pirs/src/ && rm test/Illumina_* test/ref22* test/*.out test/*.err test/EAMSS2_*
	cd ./src/stator/gcContCvgBias/src/test && rm -v test_*
	-@for P in ${PERL_LIST}; do \
	    rm `basename $${P}`; \
	done
