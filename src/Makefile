CC = gcc

OPTIMIZATIONS = -DNDEBUG -O2 -funroll-all-loops -minline-all-stringops \
				-momit-leaf-frame-pointer

CFLAGS = -W -Wformat -Wimplicit -Wreturn-type -Wall \
         -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
         -Wcast-align -Winline -Wnested-externs -Wextra \
         -std=c99

INCDIR1 = `pwd`/samtools-0.1.19
LNKDIR1 = `pwd`/samtools-0.1.19

all: asserts.h runtime.h runtime.c time.h time.c errors.h errors.c constants.h memalloc.h memalloc.c files.h files.c sequences.h sequences.c slinklist.h slinklist.c hashtable.h hashtable.c bamreport.h bamreport.c 
	cd samtools-0.1.19 && $(MAKE) && cd ..
	$(CC) $(CFLAGS) -c runtime.c
	$(CC) $(CFLAGS) -c time.c
	$(CC) $(CFLAGS) -c errors.c
	$(CC) $(CFLAGS) -c memalloc.c
	$(CC) $(CFLAGS) -c files.c	
	$(CC) $(CFLAGS) -c sequences.c
	$(CC) $(CFLAGS) -c slinklist.c
	$(CC) $(CFLAGS) -c hashtable.c
	$(CC) $(CFLAGS) -c -I$(INCDIR1) -L$(LNKDIR1) bamreport.c -lbam -lz -lm
	gcc -shared -fPIC -o libBAMreport.so *.o -L$(LNKDIR1) -lbam -lz -lm
	@echo "Add libBAMreport.so to your PYTHONPATH so that python can find it"

.PHONY : clean

clean:
	cd samtools-0.1.19 && $(MAKE) clean && cd ..
	rm *.o
	rm *.so
