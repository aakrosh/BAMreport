#ifndef BAMREPORTS_H
#define BAMREPORTS_H

void bamreports(const char* const refname,
                const char* const bamname,
                const char* const rlensname,
                const char* const qualsname,
                const char* const nucsname,
                const char* const covname,
                const char* const insname,
                const char* const gccovname,
                const int windowsize,
                const char* const statsname,
                const char* const readstarts,
                const char* const readgroups,
                const bool allinsertlengths);

#endif
