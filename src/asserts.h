#ifndef ASSERTS_H
#define ASSERTS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <execinfo.h>

#define pre  assert
#define post assert

#define forceassertwithdebug(e, f, l)    \
        if(! (e))                        \
        {                                \
            fprintf(stderr, "Assertion failed: %s file %s line %d\n", #e,f,l);\
            exit(EXIT_FAILURE);          \
        }

#define forceassert(expr) forceassertwithdebug(expr, __FILE__, __LINE__)

typedef unsigned int uint;

#endif
