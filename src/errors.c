#include "errors.h"

char* argv0;

/*print the name of the program*/
void print_argv0()
{
    pre(argv0 != NULL);

    char* p = strrchr(argv0,'/');
    fprintf(stderr,"%s: ", p ? p+1 : argv0);
}


void fatalf(const char* const fmt, ...)
{
    pre(fmt != NULL);

    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    print_argv0();
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(EXIT_FAILURE);
}

/*print the message on stderr and die*/
void fatal(const char* const msg)
{
    pre(msg != NULL);

    fatalf("%s", msg);
}

void warnf(const char* const fmt, ...)
{
    pre(fmt != NULL);

    va_list ap;
    va_start(ap, fmt);
    fflush(stdout);
    print_argv0();
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
}

/*print the message on stderr and die*/
void warn(const char* const msg)
{
    pre(msg != NULL);

    warnf("%s", msg);
}


