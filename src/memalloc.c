#include "memalloc.h"

// make a call to malloc, and panic if the allocation request fails. If enabled
// at compile time a special  memory  debugging  capabilities replaces the
// normal version of ckalloc, which aids in detecting memory overwrites and
// leaks. If the compiler flag MEM_DEBUG is set, then the debugging capabilities
// are invoked
void* ckalloc(const size_t size)
{
    void* ptr;

#if MEM_DEBUG
    // request for GUARD_SIZE bytes before and after the requested block. Write
    // NULL into that guard block. We will check that when free is called in
    // this mode.
    void* tmp;
    if((tmp = malloc(sizeof(size_t) + GUARD_SIZE + size + GUARD_SIZE)) == NULL){
        fprintf(stderr, "MEM_DEBUG : Error in allocation %zd bytes.\n", size);
        perror(NULL);
        exit(MEM_ALLOC_FAILURE);
    } 

    memset((char*)tmp + sizeof(size_t), 1, GUARD_SIZE);
    memset((char*)tmp + sizeof(size_t) + GUARD_SIZE + size, 1, GUARD_SIZE);
    *((size_t*)tmp) = size;
    
    ptr = (char*)tmp + sizeof(size_t) + GUARD_SIZE; 
#else
    if((ptr = malloc(size)) == NULL){
        fprintf(stderr, "Error in allocating %zd bytes.\n", size);
        perror(NULL);
        exit(MEM_ALLOC_FAILURE);
    }
#endif
    return ptr;
}

// fill the allocated memory with '0' before returning a pointer to it
void* ckallocz(const size_t size)
{
    void* ptr = ckalloc(size);

    if(memset(ptr, 0, size) != ptr){    
        fprintf(stderr, "Error in initializing the area\n");
        exit(MEM_SET_FAILURE);
    }

    return ptr;
}

/*reallocate the memory to the given size*/
void *ckrealloc(void * p, const size_t size)
{
    p = p ? realloc(p, size) : malloc(size);
    if (!p){
        fatal("ckrealloc failed");
    }
    return p;
}

// free the memory pointed to by this pointer. If MEM_DEBUG is set, then check
// the guards.
void ckfree(void* ptr)
{
#if MEM_DEBUG
    char* leftguard = (char*)ptr - GUARD_SIZE;
    char* szptr = leftguard - sizeof(size_t);
    size_t size = *((size_t*)szptr);  
    char* rightguard = (char*)ptr + size;

    // check the guards to make sure we dont have an error.
    int i;
    for(i = 0; i < GUARD_SIZE; i++){
        if(*(leftguard + i) != (unsigned char)1){
            fprintf(stderr, "MEM_DEBUG : Left guard tainted.\n");
            exit(MEM_ASSIGN_FAILURE);
        }
        if(*(rightguard + i) != (unsigned char)1){
            fprintf(stderr, "MEM_DEBUG : Right guard tainted.\n");
            exit(MEM_ASSIGN_FAILURE);
        }
    }
#else
    free(ptr);
#endif
}
