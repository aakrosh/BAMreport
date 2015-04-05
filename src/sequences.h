/*! Routines for reading and writing fasta sequences */

#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "asserts.h"
#include "files.h"
#include "constants.h"
#include "memalloc.h"

/*bases in a single line of sequences*/
#define FORMAT 60

/*the memory is incremented in chunks of bytes specified by this macro*/
#define CHUNK 60

/*the return value from a routine*/
typedef enum rt_status_st {ENDFILE,SUCCESS,FAILURE}rt_status;


typedef struct sequence_st
{
	uchar* header;		/*the header of the sequence*/
	uint hmax;			/*memory allocated for the header*/
	uint hlen;			/*length of the header*/

	uchar* sequence;	/*the base sequence*/
	uint smax;			/*memory allocated for the sequence*/
	uint slen;			/*length of the sequence*/
	
	FILE* fp;			/*which file do we read*/
	int rc;			    /*if set this sequence is reverse complemented*/
}sequence;

/*open the fasta file corresponding to the name of the file*/
sequence* open_fasta_sequence(const char* const file);

/*open the fasta file and return the first sequence*/
sequence* read_fasta_sequence(const char* const file);

/*pretty print this sequence*/
void format(const uchar* const sequence, const int len);

/*print the fasta sequence in a format with 60 bases in a single line, and 
 *with a fasta header*/
void print_fasta_sequence(const sequence* const sp);

/*get the next fasta sequence to the sequence sp*/
sequence* get_next_sequence(sequence* const sp); 

/*free all the resources used by this fasta sequence*/
void close_fasta_sequence(sequence* const sp);

/*clone the fasta sequence*/
sequence* copy_sequence(const sequence* const sp);

/*reverse complement the sequence in place*/
sequence* reverse_complement_sequence(sequence* const sp);

/* reverse complement this DNA sequence string in place*/
uchar* reverse_complement_string(uchar* sequence, const int length);

/*count the number of sequences in the file*/
int count_sequences(const char* const file);

/*convert the sequence of the read into a compressed format, where each
 * homopolymer stretch is replaced by a single base*/
void compress_sequence(sequence* const sp);

#endif
