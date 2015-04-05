#ifndef TIME_H
#define TIME_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "asserts.h"
#include "runtime.h"

extern time_t t0;

/* time and memory management outputs */
void timestamp(const char* const string);
void print_usage();

#endif
