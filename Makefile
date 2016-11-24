all:
	cd htslib && $(MAKE)
	cd src && $(MAKE)

.PHONY: clean install

install:
	mkdir -p bin
	cp src/BAMreport src/BAMreport.R src/bamstats bin/

clean:
	-rm -rf bin
	cd htslib && $(MAKE) clean
	cd src && $(MAKE) clean
	cd test_data && $(MAKE) clean
