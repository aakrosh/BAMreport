#ifndef MEMALLOC_H
#define MEMALLOC_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "errors.h"

#define GUARD_SIZE 8

// Routines for memory allocation
void* ckalloc(const size_t size);
void* ckallocz(const size_t size);
void* ckrealloc(void * p, const size_t size);

void ckfree(void* ptr);
#endif
