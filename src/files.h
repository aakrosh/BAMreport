#ifndef FILES_H
#define FILES_H

#include <errno.h>

#include "asserts.h"
#include "memalloc.h"

/*file handling routines*/
FILE* ckopen(const char* const name, const char* const mode);

/*drop in replacement for the getline function*/
signed long getline(char** lineptr, size_t* max, FILE* stream);

#endif
