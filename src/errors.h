#ifndef ERRORS_H
#define ERRORS_H
 
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include "asserts.h"
#include "errors.h"
  
extern char* argv0;

/*error reporting routines*/
void fatal(const char* const msg);
void fatalf(const char* const fmt, ...);

/*print the name of the program*/
void print_argv0();

#endif
