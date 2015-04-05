#ifndef SLINKLIST_H
#define SLINKLIST_H

#include <stdlib.h>

#include "asserts.h"
#include "memalloc.h"

/*abstraction for singly linked lists*/
typedef struct slinklist_st
{
	struct slinklist_st* next;
}slinklist;

#define sladdhead(plist, node) ((node)->next = *(plist), *(plist) = node)

/* add this to the tail of the list (very slow) */
void sladdtail(void* plist, const void* const node);

/*return the number of elements in the list*/
int slcount(const void* const list);

/* return the distance between the two elements in the list */
int sldistance(const void* const start, const void* const end);

/* return the ith element of the list (ith is 1-based) */
void* slelement(void* list, const uint index);

/*reverse order of  a list*/
void slreverse(void* plist);

/*remove the item from the list*/
void* slremove(void* plist, void* const node);

/*free the list and set the pointer to the list to be NULL*/
void slfreelist(void* plist);

/* return the last element of this list. Return NULL if the list is empty */
void* sllast(void* const list);

/*return and remove the first element of the list */
void* slpop(void* plist);

/* sort the linked list with qsort and a temporary array*/
void slsort(void* plist, 
		    int(*compare)(const void* const elem1, const void* const elem2));
#endif
