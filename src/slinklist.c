#include "slinklist.h"

/* add this to the tail of the list (very slow) */
void sladdtail(void* plist, const void* const node)
{
    pre(node != NULL);

    slinklist** ppt = (slinklist**) plist;
    slinklist*  n   = (slinklist*)  node;
    
    while(*ppt != NULL){
        ppt = &((*ppt)->next);
    }
    n->next = NULL;
    *ppt = n;
}


/*return the number of elements in the list*/
int slcount(const void* const list)
{
	int count = 0;
	slinklist* iter = (slinklist*) list;

	for(; iter; ++count, iter = iter->next);
	return count;
}

/* return the distance between the two elements in the list. Even if start and
 * end are the same elements, the sldistance would be 1 */
int sldistance(const void* const start, const void* const end)
{
    pre(start != NULL);

    int count = 0;
    slinklist* iter = (slinklist*)start;

    for(; iter && iter != end; ++count, iter = iter->next);
    return count+1;
}

/* return the ith element of the list (ith is 1-based). Return NULL if the index
 * is greater than the size of the list. */
void* slelement(void* list, const uint index)
{
    pre(list != NULL);

    slinklist* iter = (slinklist*) list;

    uint i = 0;
    for(; iter && i < (index - 1); i++, iter = iter->next);

    if(i == (index - 1)) return iter;
    return NULL;
}


/*reverse order of  a list*/
void slreverse(void* plist)
{
	slinklist** ppt = (slinklist**)plist;
	slinklist* newlist = NULL;
	slinklist* el = NULL;
	slinklist* next = NULL;

	next = *ppt;
	while(next != NULL){
		el = next;
		next = el->next;
		el->next = newlist;
		newlist = el;
	}
	*ppt = newlist;
}

/*remove the item from the list. */
void* slremove(void* plist, void* const node)
{
    pre(node != NULL);

	slinklist* iter = *((slinklist**) plist);
	slinklist* t = (slinklist*) node;
	slinklist* pt = (slinklist*) node;

	if(iter == t){
		*(slinklist**)plist = t->next;
		return t;
	}
	for(; iter && (iter != t); pt = iter, iter = iter->next);	
	pt->next = t->next;

    forceassert(t != NULL);
	return t;
}

/*free the list and set the pointer to the list to be NULL*/
void slfreelist(void* plist)
{
	slinklist** ppt = (slinklist**)plist;
	slinklist* next = *ppt;
	slinklist* el = NULL;

	while(next != NULL){
		el = next;
		next = el->next;
		ckfree((char*)el);
	}
	*ppt = NULL;
}

/* return the last element on the list. Return NULL if the list is empty */
void* sllast(void* const list)
{
    if(list == NULL) return NULL;

    slinklist* iter;
    for(iter = list; iter->next != NULL; iter = iter->next);
    return iter;
}

/* sort the linked list with qsort and a temporary array*/
void slsort(void* plist, 
		    int(*compare)(const void* const elem1, const void* const elem2))
{
	slinklist** pl = (slinklist**)plist;
	slinklist* list = *pl;

	int count = slcount(list);
	if(count > 1){
		slinklist* el = NULL;
		slinklist** array = NULL;
		int i = 0;

		array = ckalloc(count * sizeof(*array));
		for(el = list, i=0; el != NULL; el = el->next, i++){
			array[i] = el;
		}
		qsort(array, count, sizeof(array[0]), compare);
		list = NULL;
		for( i = 0; i < count; i++){
			array[i]->next = list;
			list = array[i];
		}
		ckfree(array);
		slreverse(&list);
		*pl = list;
	}
}

/*return and remove the first element of the list */
void* slpop(void* plist)
{
    pre(*((slinklist**) plist) != NULL);

	slinklist* node = *((slinklist**) plist);
	
	*(slinklist**)plist = node->next;
	return node;
}
