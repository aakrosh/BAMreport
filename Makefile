all:
	cd htslib && $(MAKE)
	cd libStatGen && $(MAKE)
	cd verifyBamID && $(MAKE)
	make data 
	cd src && $(MAKE)

data:
	#fetch external genotypes

.PHONY: clean install

install:
	mkdir -p bin
	cp verifyBamID/bin/verifyBamID bin/
	cp src/BAMreport src/BAMreport.R src/bamstats src/VERSION bin/ 

clean:
	-rm -rf bin
	cd htslib && $(MAKE) clean
	cd libStatGen && $(MAKE) clean
	cd verifyBamID && $(MAKE) clean
	cd src && $(MAKE) clean
	cd test_data && $(MAKE) clean
