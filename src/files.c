#include "files.h"

/*open this file and return a file pointer*/
FILE* ckopen(const char* const name, const char* const mode)
{
    pre(name != NULL);
    pre(mode != NULL);

    FILE* fp;
    if((fp = fopen(name, mode)) == NULL){
        fatalf("error in opening the file %s: %s", name, strerror(errno));
    }
    return fp;
}

/*drop in replacement for 'getline' function*/
signed long getline(char** lineptr, size_t* max, FILE* stream)
{
    int ch;
    unsigned long size = 0;
    char* ptr;
    size_t sz;

    while((ch = fgetc(stream)) != EOF){
        if(size >= *max){
            sz = size + (size >> 5) + 16;
            *max = sz;
            if((ptr = ckrealloc(*lineptr, *max)) == NULL){
                return -1;
            }
            *lineptr = ptr;
        }
        (*lineptr)[size++] = ch;
        if(ch == '\n'){
            break;
        }
    }

    if(size != 0){
        if(size >= *max){
            sz = size + (size >> 5) + 16;
            *max = sz;
            if((ptr = ckrealloc(*lineptr, *max)) == NULL){
                return -1;
            }
            *lineptr = ptr;
        }
        (*lineptr)[size] = '\0';
    }
    if(0 == size || ch == EOF){
        return -1;
    }

   return size;
}

