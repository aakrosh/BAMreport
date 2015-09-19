CC = gcc

OPTIMIZATIONS = -DNDEBUG -O2 -minline-all-stringops -momit-leaf-frame-pointer

CFLAGS = -W -Wformat -Wimplicit -Wreturn-type -Wall \
         -Wunused-variable -Wunused-parameter -Wreturn-type -Wswitch \
         -Winline -Wnested-externs -Wextra \
         -std=c99

INCDIR1 = $(PWD)/htslib/htslib
LNKDIR1 = $(PWD)/htslib

ifneq ($(shell uname -s),Darwin)
	OPTIMIZATIONS += -funroll-all-loops
endif 

DEBUG ?= 0

ifeq ($(DEBUG), 1)
	CFLAGS += -g
else
	CFLAGS += $(OPTIMIZATIONS)
endif

all: bamstats

bamstats: asserts.h runtime.h runtime.c time.h time.c errors.h errors.c constants.h memalloc.h memalloc.c files.h files.c sequences.h sequences.c slinklist.h slinklist.c hashtable.h hashtable.c bamstats.c
	cd htslib && $(MAKE)
	$(CC) $(CFLAGS) -c runtime.c
	$(CC) $(CFLAGS) -c time.c
	$(CC) $(CFLAGS) -c errors.c
	$(CC) $(CFLAGS) -c memalloc.c
	$(CC) $(CFLAGS) -c files.c      
	$(CC) $(CFLAGS) -c sequences.c
	$(CC) $(CFLAGS) -c slinklist.c
	$(CC) $(CFLAGS) -c hashtable.c
	$(CC) $(CFLAGS) -I$(INCDIR1) -L$(LNKDIR1) -o bamstats *.o bamstats.c htslib/libhts.a -lz -lm -lpthread 

.PHONY : clean

clean:
	cd htslib && $(MAKE) clean && cd ..
	-rm *.o
	-rm *.so
	-rm bamstats

