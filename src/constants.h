#ifndef CONSTANTS_H
#define CONSTANTS_H

// bool constants
typedef enum bool_en{
    TRUE  = 1,
    FALSE = 0
}bool;

typedef unsigned char uchar;

#undef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))

// exit codes
#define ASSERT_FAILURE 2
#define MEM_ALLOC_FAILURE 3
#define MEM_ASSIGN_FAILURE 4
#define MEM_SET_FAILURE 5

#endif
