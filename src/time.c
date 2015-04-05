#include "time.h"

time_t t0;

/* time and memory management outputs */
void timestamp(const char* const string)
{
    pre(string != NULL);

    fprintf(stderr, "%s", string);
    time_t t1 = time(0);
    fprintf(stderr, "%3.2f sec. elapsed\n", difftime(t1,t0));
    fprintf(stderr, "--------------------------------------------\n");
}

void print_usage()
{
    struct rusage usage;
    if(getrusage(RUSAGE_SELF, &usage) == 0){
        fprintf(stderr, "Program stats:\n");
        fprintf(stderr, "User time used  : %ld\n", usage.ru_utime.tv_sec);    
        fprintf(stderr, "System time used: %ld\n", usage.ru_stime.tv_sec); 
        fprintf(stderr, "Maximum resident set size: %ld bytes\n", 
        usage.ru_maxrss);   
        fprintf(stderr, "--------------------------------------------\n");
    }
}
