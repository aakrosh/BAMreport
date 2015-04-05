#include "hashtable.h"

/*allocate a new hash table*/
hashtable* new_hashtable(const int po2size)
{
	int max_size = po2size;
	if(po2size > 24){
		fprintf(stderr,"hash po2size should not exceed 24, using 24\n");
		max_size = 24;
	}

	hashtable* hashtable = ckallocz(sizeof(struct hashtable_st));
	hashtable->po2size = max_size;
	hashtable->size = (1 << hashtable->po2size);
	hashtable->mask = hashtable->size - 1;

	hashtable->bins  = ckallocz(hashtable->size*sizeof(struct bin_st*));

	return hashtable;
}

#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
  || defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((const uint8_t *)(d))[1] << UINT32_C(8))\
                      +((const uint8_t *)(d))[0])
#endif

static uint32_t superfasthash (const char* data, int len) {
uint32_t hash = len, tmp;
int rem;

    if (len <= 0 || data == NULL) return 0;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
                hash ^= hash << 16;
                hash ^= data[sizeof (uint16_t)] << 18;
                hash += hash >> 11;
                break;
        case 2: hash += get16bits (data);
                hash ^= hash << 11;
                hash += hash >> 17;
                break;
        case 1: hash += *data;
                hash ^= hash << 10;
                hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}


/*add the following substring of length len to the hashtable. Return the 
 * bin corresponding to it*/
bin* add_hashtable(hashtable* const hash,  /*the hashtable*/
				   const char* const name, /*the string*/
				   const int length, 	   /*use 'length' characters of name*/
				   void* val)				/*the value*/
{
    pre(name != NULL);

	uint32_t index =  superfasthash(name, length);
	index = index & hash->mask;
	
	char* hname = ckalloc(length+1);
	memcpy(hname, name, length);
	hname[length] = '\0';

	bin* bin = ckallocz(sizeof(struct bin_st));
	bin->name = hname;
	bin->val = val;
	if(hash->bins[index]){
		hash->collisions++;
	}
	bin->next = hash->bins[index];
	hash->bins[index] = bin;
	hash->elcount++;
	
    post(bin != NULL);
	return bin;
}

/*look up the hash table to see an entry corresponding to the string of length
 * len exists. Return the bin if it does exist*/
void* lookup_hashtable(hashtable* const hash, 
				       const char* const name, 
				       const int len)
{
    pre(name != NULL);

	uint32_t index = superfasthash(name, len);
	index = index & hash->mask;
	
	bin* iter = NULL;
	bin* bin = NULL;
	if(hash->bins[index] != NULL){
		for(iter = hash->bins[index];iter; iter = iter->next){
			if(strncmp(iter->name, name, len) == 0){
				bin = iter;
			}
		}
	}
	return bin;
}

/*look in the hashtable and return the value associated with the string of
 * length len. Returns NULL is the bin doesnt exist*/	
void* find_value_hashtable(hashtable* const hash,
						   const char* const name,
						   const int len)
{
    pre(name != NULL);

	bin* bin = lookup_hashtable(hash, name, len);
	if(bin == NULL){
		return NULL;
	}
	return bin->val;
}

/* find a value in the hashtable or die if the key does not exist in the hash*/
void* must_find_hashtable(hashtable* const hash,
						  const char* const name,
						  const int len)
{
    pre(name != NULL);

	bin* bin = lookup_hashtable(hash, name, len);
	if(bin == NULL){
		fatalf("did not find %s in the hash", name);
	}
	return bin->val;
}

/*remove this entry from the hashtable and return the value associated with it*/
void* remove_hashtable_entry(hashtable* const hash,
							 const char* const name,
							 const int len)
{
    pre(name != NULL);

	uint32_t index = superfasthash(name, len);
	index = index & hash->mask;
	
	bin* iter = NULL;
	bin* bin = NULL;
	if(hash->bins[index] != NULL){
		for(iter = hash->bins[index];iter; iter = iter->next){
			if(strncmp(iter->name, name, len) == 0){
				bin = iter;
				break;
			}
		}
	}

	void* val = NULL;

	if(bin){
		val = bin->val;
		slremove(&hash->bins[index], bin);
		ckfree(bin->name);
		ckfree(bin);
		hash->elcount--;
	}
	return val;
}

/* apply this function to every non-null bin in the hashtable */
void func_hashtable(hashtable* const hash, 
					void (*func)(bin* bin))
{
    pre(func != NULL);

	int i;
	bin* iter;
	bin* next;
	
	for(i = 0; i < hash->size; i++){
		iter = hash->bins[i];
		while(iter){
			next = iter->next;
			func(iter);
			iter = next;
		}
	}
}

static void free_hash(bin* bin)
{
	ckfree(bin->name);
	ckfree((char*)bin);
}

static void comp_free_hash(bin* bin)
{
	ckfree(bin->name);
	ckfree(bin->val);
	ckfree((char*)bin);
}

/*free all the resources allocated to the hash table*/
void free_hashtable(hashtable** const phash)
{
	if(*phash != NULL){
		func_hashtable(*phash, free_hash);

		hashtable* hash = *phash;
		ckfree(hash->bins);
		ckfree(hash);
		*phash = NULL;
	}
}

/*free all the resources allocated to the hash table. Also for each bin, free
 * the value in the bin.*/
void free_hashtable_completely(hashtable** const phash)
{	
	if(*phash != NULL){
		func_hashtable(*phash, comp_free_hash);

		hashtable* hash = *phash;
		ckfree(hash->bins);
		ckfree(hash);
		*phash = NULL;
	}
}

static void print_names(bin* bin)
{
    pre(bin != NULL);
	printf("%s\n", bin->name);
}

/* for debugging only. print all the names in the hashtable */
void print_hashtable(hashtable* const hash)
{
	func_hashtable(hash, print_names);
}

