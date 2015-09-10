CC = gcc

OPTIMIZATIONS = -DNDEBUG -O2 -funroll-all-loops -minline-all-stringops \
				-momit-leaf-frame-pointer

CFLAGS = -W -Wformat -Wimplicit -Wreturn-type -Wall \
         -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
         -Wcast-align -Winline -Wnested-externs -Wextra \
         -std=c99 -g

INCDIR1 = `pwd`/samtools-0.1.19
LNKDIR1 = `pwd`/samtools-0.1.19

all: asserts.h runtime.h runtime.c time.h time.c errors.h errors.c constants.h memalloc.h memalloc.c files.h files.c sequences.h sequences.c slinklist.h slinklist.c hashtable.h hashtable.c bamstats.c 
	cd samtools-0.1.19 && $(MAKE) && cd ..
	$(CC) $(CFLAGS) -c runtime.c
	$(CC) $(CFLAGS) -c time.c
	$(CC) $(CFLAGS) -c errors.c
	$(CC) $(CFLAGS) -c memalloc.c
	$(CC) $(CFLAGS) -c files.c	
	$(CC) $(CFLAGS) -c sequences.c
	$(CC) $(CFLAGS) -c slinklist.c
	$(CC) $(CFLAGS) -c hashtable.c
	$(CC) $(CFLAGS) -I$(INCDIR1) -L$(LNKDIR1) -o bamstats *.o bamstats.c -lbam -lz -lm -lpthread

.PHONY : clean

clean:
	cd samtools-0.1.19 && $(MAKE) clean && cd ..
	-rm *.o
	-rm *.so
	-rm bamstats
