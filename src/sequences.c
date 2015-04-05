/*! Routines for reading and writing fasta file */

#include "sequences.h"

/*all valid bases in the ASCII world should be 1*/
static const char bases[] = {
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,1,1,1,1,0,
0,1,1,0,0,1,0,1,1,0,
0,0,1,1,1,0,1,1,0,1,
0,0,0,0,0,0,0,1,1,1,
1,0,0,1,1,0,0,1,0,1,
1,0,0,0,1,1,1,0,1,1,
0,1,0,0,0,0,0,0,0
};

static const unsigned char dna_complement[] =
  "                                                                "
  " TVGH  CD  M KN   YSA BWXR       tvgh  cd  m kn   ysa bwxr      "
  "                                                                "
  "                                                                ";


/*open the fasta file corresponding to the name of the file*/
sequence* open_fasta_sequence(const char* const file)
{   
    pre(file != NULL);

	sequence* sp = ckallocz(sizeof(struct sequence_st));
	FILE* fp = ckopen(file,"r");

	sp->fp = fp;
	sp->rc = FALSE;
	return sp;
}

void append(uchar** parray,     /*the pointer to the array to be used*/
            uint* const pmax,   /*how many bytes have been allocated*/
            uint* const plen,   /*how many bytes have been used*/
            const int ch)           /*the byte to be added*/
{
    uint max = *pmax;
    uint len = *plen;

    /*do we need more memory*/
    if(len >= max){
        *parray = ckrealloc(*parray, max + CHUNK);
        *pmax = *pmax + CHUNK;
    }
    (*parray)[len] = ch;
    *plen = *plen + 1;
}


/*read the sequence*/
static rt_status read_sequence(sequence* const sp)
{
    pre(sp != NULL);

	/*is it the end of the file?*/
	if(feof(sp->fp) || ferror(sp->fp)){
		return ENDFILE;
	}
	
	int ch;

	/*skip all the whitespace characters*/
	do{
		ch = fgetc(sp->fp);
	}while(ch == ' ' || ch == '\t');

	/*check and add the fasta header*/
	assert(ch == '>');
	ch = fgetc(sp->fp);
	do{
		append(&sp->header, &sp->hmax, &sp->hlen, ch);
		ch = fgetc(sp->fp);
	}while(ch != '\n');
	append(&sp->header, &sp->hmax, &sp->hlen, 0);
	sp->hlen--;

	/*add the bases*/
	do{
		ch = fgetc(sp->fp);
		if(1 == bases[ch]){
			append(&sp->sequence, &sp->smax, &sp->slen, ch);
		}
	}while(ch != '>' && ch != EOF);
	append(&sp->sequence, &sp->smax, &sp->slen, 0);
	sp->slen--;

	if(ch != EOF && ungetc(ch, sp->fp) == EOF){
		fatalf("error in pushing back the >");
	}
	
	return SUCCESS;
}

/*open the fasta file and return the first sequence*/
sequence* read_fasta_sequence(const char* const file)
{
    pre(file != NULL);

	/*open the fasta file*/
	sequence* sp = open_fasta_sequence(file);

	/*read the sequence*/
	if(read_sequence(sp) == ENDFILE){
		return NULL;
	}

    post(sp != NULL);
	return sp;
}

/*free the resources with this sequence*/
void close_fasta_sequence(sequence* sp)
{
    if(sp != NULL){
        ckfree(sp->header); 
        ckfree(sp->sequence);
        if(sp->fp != NULL && sp->fp != stdin){  
            fclose(sp->fp);
        } 
        ckfree(*(&sp));
        *(&sp) = NULL;
    }
}

/*get the next fasta sequence to the sequence sp*/
sequence* get_next_sequence(sequence* const sp)
{
    pre(sp != NULL);

    sequence* psp = sp;
	sp->slen = 0;
	sp->hlen = 0;
	if(read_sequence(sp) == ENDFILE){
        close_fasta_sequence(psp);
		return NULL;
	}

	return sp;
}

/*pretty print this sequence*/
void format(const uchar* const sequence, const int len)
{
    pre(sequence != NULL);

	int i,j;
	for(i = 0, j = 0; i < len; i++){
		if(j != 0  && (j % 60)==0){
			printf("\n");
		}
		if(sequence[i] != '-'){
			printf("%c", sequence[i]);
			j++;
		}
	}
	printf("\n");
}

/*print the fasta sequence in a format with 60 bases in a single line, and 
 *with a fasta header*/
void print_fasta_sequence(const sequence* const sp)
{
    pre(sp != NULL);

	printf(">%s\n", sp->header);
	format(sp->sequence, sp->slen);
}

static char* copy_string(const char* const str)
{
    pre(str != NULL);

    char* copy = ckalloc(strlen(str)+1);
    strcpy(copy,str);
    return copy;
}

/*clone the fasta sequence*/
sequence* copy_sequence(const sequence* const sp)
{
    pre(sp != NULL);

	sequence* sequence = ckallocz(sizeof(struct sequence_st));

	sequence->hmax = sequence->hlen = strlen((char*)sp->header);
	sequence->header = (uchar*)copy_string((char*)sp->header);
	
	sequence->smax = sequence->slen = strlen((char*)sp->sequence);
	sequence->sequence = (uchar*)copy_string((char*)sp->sequence);
	return sequence;
}

/*reverse complement the string in place*/
uchar* reverse_complement_string(uchar* sequence, const int length)
{
    pre(sequence != NULL);

	uchar* s = sequence;
	uchar* p = s + length - 1;

	while(s <= p){
		uchar c;

		c = dna_complement[*s];
		*s = dna_complement[*p];
		*p = c;
		++s; --p;
	}
	return sequence;
}

/*reverse complement the sequence in place*/
sequence* reverse_complement_sequence(sequence* const sp)
{
    pre(sp != NULL);

	sp->sequence = reverse_complement_string(sp->sequence, sp->slen);	
	sp->rc = sp->rc == TRUE ? FALSE : TRUE;
	return sp;
}

/*count the number of sequences in the file*/
int count_sequences(const char* const file)
{
    pre(file != NULL);

	sequence* sp;
	int count = 0;

	if((sp = read_fasta_sequence(file)) == NULL){
		fatalf("error in reading the fasta reads from %s", file);
	}

	while(sp){
		count++;
		if(!get_next_sequence(sp)){
			break;
		}
	}
	close_fasta_sequence(sp);

	return count;
}

/*convert the sequence of the read into a compressed format, where each
 * homopolymer stretch is replaced by a single base*/
void compress_sequence(sequence* const sp)
{
    pre(sp != NULL);

	char base;
	uint i = 0, j = 0;
	
	while(i < sp->slen){
		base = sp->sequence[i];
		while((i+1) < sp->slen && sp->sequence[i+1] == base){
			i++;
		}
		i++;
		sp->sequence[j++] = base; 
	}
	sp->sequence[j] = '\0';
	sp->slen = j;
}
