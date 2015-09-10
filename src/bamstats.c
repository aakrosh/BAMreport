#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <inttypes.h>

#include "bam.h"

#include "errors.h"
#include "time.h"
#include "memalloc.h"
#include "files.h"
#include "sequences.h"
#include "hashtable.h"

// This program is useful to generate information from a BAM file. The following
// information is provided as output:
// Following statistics about the alignment:
//
// QC-passed reads : 12677423
// QC-failed reads : 0
// Number of reads paired in sequencing: 0
// Number of single-fragments: 12677423
// 
// Duplicate reads : 3973834 (31.35%)
// Number of reads aligned: 6134567 (48.39%)
// Number of reads mapped: 4619967 (36.44%)
// 
// Number of read1's : 0
// Number of read2's : 0
// Number of reads properly paired: 0(nan)
// Number of reads with itself and mate aligned: 0
//
// QC-passed bases : 
// QC-failed bases : 0
// Number of bases paired in sequencing:
// Number of single-fragment bases     : 
// 
// Duplicate bases : 
// Number of bases aligned : 
// Number of bases mapped  : 
//
// Besides this we generate (if the user requests)
// a) quality value distribution with position on read.
// b) nucleotide composition for the reads. 
// c) coverage distribution of the reference  sequence.
// d) insert distance distribution.
// e) GC based coverage.

// print debug information?
int debug_flag = FALSE;

// the maximum quality value allowed
const uint maxqualityvalue = 60;

bool nowarnings;

typedef struct bamstats_t{
    // generation statistics
    uint64_t r_qcpassed;
    uint64_t r_qcfailed;
    uint64_t r_paired;
    uint64_t r_singletons;
    uint64_t b_qcpassed;
    uint64_t b_qcfailed;
    uint64_t b_paired;
    uint64_t b_singletons;
    
    // alignment statistics
    uint64_t r_duplicates;
    uint64_t r_aligned;
    uint64_t r_mapped;
    uint64_t b_duplicates;
    uint64_t b_aligned;
    uint64_t b_mapped;

    uint64_t read1s;
    uint64_t read2s;
    uint64_t properlypaired;
    uint64_t bothaligned;
}bamstats;

typedef struct chrcoverage_t
{
    uint64_t length;
    uchar* coverage;
    uchar* sequence;
} chrcoverage;

// return the number of bases in the read that are used in the alignment
static int lookatcigar(const bam1_t* const alignment)
{
    int numused = 0;

    uint32_t* cigar = bam1_cigar(alignment);

    uint i;
    for(i = 0; i < alignment->core.n_cigar; i++){
        switch(cigar[i] & 0xf){
            case BAM_CMATCH :
                numused += ((cigar[i] & 0xfffffff0) >> 4); break;
            case BAM_CINS:
                numused += ((cigar[i] & 0xfffffff0) >> 4); break;
            case BAM_CDEL:
                break;
            case BAM_CREF_SKIP:
                fatalf("Have not handled REF_SKIP in CIGAR");
            case BAM_CSOFT_CLIP:
                break;
            case BAM_CHARD_CLIP:
                break;
            case BAM_CPAD:
                break;
            default :
                fatalf("Unhandled CIGAR operation\n");
        }
    }

    return numused;
}

// print the statistics for this BAM file
static void print_stats(const bamstats* stats, FILE* fp)
{
    fprintf(fp, "Reads Summary\n");
    fprintf(fp,"-------------\n");
    fprintf(fp,"QC-passed reads : %"PRIu64"\n", stats->r_qcpassed);
    fprintf(fp,"QC-failed reads : %"PRIu64"\n", stats->r_qcfailed);
    fprintf(fp,"Number of single-fragments           : %"PRIu64"\n", stats->r_singletons);
    fprintf(fp,"Number of reads paired in sequencing : %"PRIu64"\n", stats->r_paired);
    fprintf(fp,"\n");
    fprintf(fp,"Duplicate reads         : %"PRIu64" (%2.2f%%)\n", 
            stats->r_duplicates, 
            stats->r_duplicates * 100.0 / stats->r_qcpassed);
    fprintf(fp,"Number of reads aligned : %"PRIu64" (%2.2f%%)\n", 
            stats->r_aligned,
            stats->r_aligned * 100.0 / stats->r_qcpassed);
    fprintf(fp,"Number of reads mapped  : %"PRIu64" (%2.2f%%)\n", 
            stats->r_mapped,
            stats->r_mapped * 100.0 / stats->r_qcpassed);
    fprintf(fp,"\n");
    fprintf(fp,"Number of read1's : %"PRIu64"\n", stats->read1s);
    fprintf(fp,"Number of read2's : %"PRIu64"\n", stats->read2s);
    fprintf(fp,"Number of reads properly paired              : %"PRIu64"\n", 
            stats->properlypaired);
    fprintf(fp,"Number of reads with itself and mate aligned : %"PRIu64"\n", 
            stats->bothaligned);
    fprintf(fp,"\n");
    fprintf(fp,"Bases Summary\n");
    fprintf(fp,"-------------\n");
    fprintf(fp,"QC-passed bases : %"PRIu64"\n", stats->b_qcpassed);
    fprintf(fp,"QC-failed bases : %"PRIu64"\n", stats->b_qcfailed);
    fprintf(fp,"Number of bases in single-fragments  : %"PRIu64"\n", stats->b_singletons);
    fprintf(fp,"Number of bases paired in sequencing : %"PRIu64"\n", stats->b_paired);
    fprintf(fp,"\n");
    fprintf(fp,"Duplicate bases         : %"PRIu64" (%2.2f%%)\n", 
            stats->b_duplicates, 
            stats->b_duplicates * 100.0 / stats->b_qcpassed);
    fprintf(fp,"Number of bases aligned : %"PRIu64" (%2.2f%%)\n", 
            stats->b_aligned,
            stats->b_aligned * 100.0 / stats->b_qcpassed);
    fprintf(fp,"Number of bases mapped  : %"PRIu64" (%2.2f%%)\n", 
            stats->b_mapped,
            stats->b_mapped * 100.0 / stats->b_qcpassed);
    fprintf(fp,"\n");
}

// print the read length distribution 
static void print_rlength(const char* const rlensname,
              const uint64_t* const rlens1,
              const uint64_t* const rlens2,
              const uint maxrlen)
{
    uint i;
    if(rlensname != NULL){
        FILE* fp = ckopen(rlensname, "w");

        for(i = 0; i < maxrlen; i++){
            fprintf(fp, "1 %d %"PRIu64"\n", i + 1, rlens1[i]);
        }
        for(i = 0; i < maxrlen; i++){
            fprintf(fp, "2 %d %"PRIu64"\n", i + 1, rlens2[i]);
        }

        fclose(fp);
    }
}

static void add_qualitycounts(const bam1_t* const alignment,
                              uint64_t** const qualvalues)
{
    uchar* qual = bam1_qual(alignment);
//    if((*qual) == '*'){
//        if(nowarnings == FALSE){
//            fprintf(stderr,"The quality values are not stored in this BAM for read %s\n", bam1_qname(alignment));
//        }
//        return;
//    }

    int i;
    for(i = 0;  i < alignment->core.l_qseq; i++){
        qualvalues[i][qual[i]>=maxqualityvalue ? maxqualityvalue-1 : qual[i]] += 1;
    }
}

void summary(const uint64_t* const frequencytable,    
             int* const pwmin, 
             int* const pbmin,
             int* const pmedian,
             int* const pbmax, 
             int* const pwmax)
{
    // find the number of entries 
    uint64_t numentries = 0;
    for(uint i = 0; i < maxqualityvalue; i++){
        numentries += frequencytable[i];
    }

    // the index of the 25th quartile and the 75th quartile 
    int index25 = numentries/4;
    int index50 = numentries/2;
    int index75 = 3*numentries/4;

    *pwmin = 49;
    *pwmax = 0;
    *pbmin   = -1;
    *pmedian = -1;
    *pbmax   = -1;

    int i = 0;
    int count  = 0;
    for(i = 0; i < (int)maxqualityvalue; i++){
        if(frequencytable[i] > 0){
            count += frequencytable[i];
            if(count >= index25 && *pbmin == -1){
                *pbmin = i;    
            }
            if(count >= index50 && *pmedian == -1){
                *pmedian = i;
            }
            if(count >= index75 && *pbmax == -1){
                *pbmax = i;
            }

            if(i < *pwmin){
                *pwmin = i;    
            }
            if(i > *pwmax){
                *pwmax = i;
            }
        }
    }
}



// print the quality value distribution
static void print_readqualdist(const char* const qualsname,     
                               uint64_t** const qual1values, 
                               uint64_t** const qual2values, 
                               const uint maxreadlength)
{
    uint i;
    if(qualsname != NULL){
        FILE* fp = ckopen(qualsname, "w");

        int whisker_min, box_min, median, box_high, whisker_high;
        for(i = 0; i < maxreadlength; i++){
            summary(qual1values[i],
                    &whisker_min, &box_min, &median, &box_high, &whisker_high);
            fprintf(fp, "1 %d %d %d %d %d %d\n",
                    i+1, whisker_min, box_min, median, box_high, whisker_high);
        }
        for(i = 0; i < maxreadlength; i++){
            summary(qual2values[i],
                    &whisker_min, &box_min, &median, &box_high, &whisker_high);
            fprintf(fp, "2 %d %d %d %d %d %d\n",
                    i+1, whisker_min, box_min, median, box_high, whisker_high);
        }
        fclose(fp);
    }
}

static void print_nuccomp(const char* const nucsname, 
                          uint64_t** const nuc1counts, 
                          uint64_t** const nuc2counts, 
                          const uint maxreadlength)
{
    uint i;
    if(nucsname != NULL){
        FILE* fp = ckopen(nucsname, "w");
    
        for(i = 0; i < maxreadlength; i++){
            fprintf(fp, "1 %u %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n",
                    i + 1,
                    nuc1counts[i][0],
                    nuc1counts[i][1],
                    nuc1counts[i][2],
                    nuc1counts[i][3],
                    nuc1counts[i][4],
                    nuc1counts[i][0] +
                    nuc1counts[i][1] + 
                    nuc1counts[i][2] +
                    nuc1counts[i][3] + 
                    nuc1counts[i][4]);
        }
        for(i = 0; i < maxreadlength; i++){
            fprintf(fp, "2 %d %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64" %"PRIu64"\n",
                    i + 1,
                    nuc2counts[i][0],
                    nuc2counts[i][1],
                    nuc2counts[i][2],
                    nuc2counts[i][3],
                    nuc2counts[i][4],
                    nuc2counts[i][0] +
                    nuc2counts[i][1] + 
                    nuc2counts[i][2] +
                    nuc2counts[i][3] + 
                    nuc2counts[i][4]);
        }
        fclose(fp);
    }
   
}

// add the nucleotide counts from this read to the given counters
static void add_nucleotidecounts(const bam1_t* const alignment,
                                 uint64_t** const nuccounts)
{
    uchar* seq = bam1_seq(alignment);
    if(*seq == '*'){
        fatalf("The read sequences is not stored in this BAM for this read %s\n", bam1_qname(alignment));
    }
    int i;
    for(i = 0; i < alignment->core.l_qseq; i++){
        switch(bam1_seqi(seq, i)){
            case 1 : nuccounts[i][0] += 1; break;
            case 2 : nuccounts[i][1] += 1; break;
            case 4 : nuccounts[i][2] += 1; break;
            case 8 : nuccounts[i][3] += 1; break;
            case 15: nuccounts[i][4] += 1; break;
            default: 
                fatalf("Incorrect base call.\n");
        }
    }
}

static void add_clippingcounts(const bam1_t* const alignment, 
                               uint64_t* const clipfrequency)
{
    uint32_t* cigar = bam1_cigar(alignment);
    bool is_reverse = (alignment->core.flag & 0x10) == 0x10 ? TRUE : FALSE;
    uint rlen = alignment->core.l_qseq;

    // only the 5' soft-clips matter for this.
    uint indx = 0;
    if (is_reverse) indx = alignment->core.n_cigar-1;
        
    if (((cigar[indx] & 0xf) == BAM_CSOFT_CLIP) || 
        ((cigar[indx] & 0xf) == BAM_CHARD_CLIP)) {
        clipfrequency[((cigar[indx] & 0xfffffff0) >> 4)] += 1;
    }
}

// print the insert length distribution for the reads
static void print_insertlendist(const char* const insname, 
                                const uint64_t* const insertlenfrequency, 
                                const uint maxinsertsize)
{
    if(insname != NULL){
        FILE* fp = ckopen(insname, "w");

        uint i;
        for(i = 0; i < maxinsertsize; i++){
            fprintf(fp, "%u %"PRIu64"\n", i + 1, insertlenfrequency[i]);    
        }

        fclose(fp);
    }
}

// print the clipping frequency for the sequences
static void print_clippingcounts(const char* const clipname,
                                 const uint64_t* const clip1frequency, 
                                 const uint64_t* const clip2frequency,
                                 const uint64_t maxrlen)
{
    uint i;
    if(clipname != NULL) {
        FILE* fp = ckopen(clipname, "w");

        for(i = 0; i < maxrlen; i++){
            fprintf(fp, "1 %d %"PRIu64"\n", i + 1, clip1frequency[i]);
        }
        for(i = 0; i < maxrlen; i++){
            fprintf(fp, "2 %d %"PRIu64"\n", i + 1, clip2frequency[i]);
        }

        fclose(fp);
    }
}

// read the reference sequence. Return a hashtable where the key is the name of
// the chromosome and the value is a coverage array (initialized to 0) with a
// size equal to the size of the chromosome.
static hashtable* read_reference(const char* const refname)
{
    hashtable* reference = new_hashtable(12);

    sequence* sp = read_fasta_sequence(refname);
    
    while(sp != NULL){
        // allocate a coverage array for the sequence
        chrcoverage* cov = ckalloc(sizeof(chrcoverage));
        cov->length   = strlen((char*)sp->sequence);
        cov->coverage = ckallocz(strlen((char*)sp->sequence));
        cov->sequence = ckallocz(strlen((char*)sp->sequence)+1);
        cov->sequence = (uchar*)strcpy((char*)cov->sequence, (char*)sp->sequence);

        // if the name of the sequence has more than one tokens, just use the
        // first token in the name
        int i = 0;
        while((sp->header[i] != '\n') && 
              (sp->header[i] != 0)    && 
              (sp->header[i] != '\t') && 
              (sp->header[i] != 32)) i++;
        sp->header[i] = 0;

        add_hashtable(reference, (char*)sp->header, strlen((char*)sp->header), cov);
        sp = get_next_sequence(sp);
    } 

    return reference;
}

// print the coverage distribution of the reference sequence
static void print_coverage_distribution(const char* const covname, 
                                        const hashtable* const reference)
{
    if(covname == NULL) return;

    uint* coverage = ckallocz(251 * sizeof(uint));

    int i;
    bin* iter;
    bin* next;
    for(i = 0; i < reference->size; i++){
        iter = reference->bins[i];
        while(iter){
            next = iter->next;
         
            uint j;       
            chrcoverage* tmp = iter->val;
            for(j = 0; j < tmp->length; j++){
                assert(tmp->coverage[j] <= 250);
                coverage[tmp->coverage[j]] += 1;
            }
            
            iter = next;
        }
    }

    FILE* fp = ckopen(covname, "w");
    
    for(i = 0; i < 251; i++){
        fprintf(fp, "%d %u\n", i, coverage[i]);
    }

    fclose(fp);
    ckfree(coverage);
}

// print the information
// chromosome start stop GC(in fraction) avgcoverage num_nonN_bases
static void  print_gccoverage(const char* const gccovname, 
                              const hashtable* const reference,
                              const uint windowsize)
{
    if(gccovname == NULL) return;

    FILE* fp = ckopen(gccovname, "w");

    int i;
    uint j, k;
    int a,c,g,t,n,count;
    bin* iter;
    bin* next;
    for(i = 0; i < reference->size; i++){
        iter = reference->bins[i];
        while(iter){
            next = iter->next;
                
            chrcoverage* tmp = iter->val;
            
            for(j = 0; j < tmp->length; j += windowsize){
                a = c = g = t = n = count = 0;
                for(k = j; k < (j + windowsize) && k < tmp->length; k++){
                    switch(tmp->sequence[k]){
                        case 'A': a += 1; count += tmp->coverage[k]; break;
                        case 'a': a += 1; count += tmp->coverage[k]; break;
                        case 'C': c += 1; count += tmp->coverage[k]; break;
                        case 'c': c += 1; count += tmp->coverage[k]; break;
                        case 'G': g += 1; count += tmp->coverage[k]; break;
                        case 'g': g += 1; count += tmp->coverage[k]; break;
                        case 'T': t += 1; count += tmp->coverage[k]; break;
                        case 't': t += 1; count += tmp->coverage[k]; break;    
                        case 'N': n += 1; break;
                        case 'n': n += 1; break;     
                        default : break; 
                    }
                }

                if((a+c+g+t) != 0){
                    fprintf(fp, "%s %d %"PRIu64" %0.3f %2.2f %d\n", 
                    iter->name, j, 
                    (j+windowsize) < tmp->length ? j+windowsize : tmp->length,
                    (g+c)*1.0/(a+c+g+t), count*1.0/(a+c+g+t), a+c+g+t);
                }
            }
            
            iter = next;
        }
    }
    fclose(fp);
}


void calcbamstats(const char* const refname,
                  const char* const bamname,
                  const char* const rlensname,
                  const char* const qualsname,
                  const char* const nucsname,
                  const char* const covname,
                  const char* const insname,
                  const char* const gccovname,
                  const char* const clipname,
                  const int windowsize,
                  const char* const statsname,
                  const char* const readstarts,
                  const char* const readgroups,
                  const bool allinsertlengths)
{
    // lets read the reference sequence. The reference sequence will be a
    // hashtable where the key is the name of the chromosome and the value will
    // be a char array of size equal to  the length of the chromosome sequence.
    hashtable* reference = NULL;
    if((covname != NULL) || (gccovname != NULL)){
        reference = read_reference(refname);
        if(nowarnings == FALSE){
            //fprintf(stderr, "Warning : reference sequence in memory (expensive)\n");    
        }
        //timestamp("Read the reference sequences\n");
    }

    int i; 

    // open the BAM file and read and process one record at a time.
    bamFile fp = bam_open(bamname, "r");

    // read the header from the BAM file
    bam_header_t* hin = bam_header_read(fp);

    // structure to hold one alignment record
    bam1_t* alignment = ckallocz(sizeof(bam1_t));

    // the statistics of the alignment will be stored in this.
    bamstats* stats = ckallocz(sizeof(bamstats));

    // the maximum read length seen so far
    int maxreadlength = 1;

    // store the read length frequencies here
    uint64_t* readlen1frequency = ckallocz(maxreadlength * sizeof(uint64_t));
    uint64_t* readlen2frequency = ckallocz(maxreadlength * sizeof(uint64_t));

    // store the quality value distribution here
    uint64_t** qual1values = ckallocz(maxreadlength * sizeof(uint64_t*));
    uint64_t** qual2values = ckallocz(maxreadlength * sizeof(uint64_t*));
    for(i = 0; i < maxreadlength; i++){
        qual1values[i] = ckallocz(maxqualityvalue * sizeof(uint64_t));
        qual2values[i] = ckallocz(maxqualityvalue * sizeof(uint64_t));
    }

    // store the nucleotide distribution here
    uint64_t** nuc1counts = ckallocz(maxreadlength * sizeof(uint64_t*));
    uint64_t** nuc2counts = ckallocz(maxreadlength * sizeof(uint64_t*));
    for(i = 0; i < maxreadlength; i++){
        nuc1counts[i] = ckallocz(5 * sizeof(uint64_t));
        nuc2counts[i] = ckallocz(5 * sizeof(uint64_t));
    }

    // store the insert length distribution here
    int maxinsertsize = 1;
    uint64_t* insertlenfrequency = ckallocz(maxinsertsize * sizeof(uint64_t));
    //timestamp("Allocated the required memory\n");

    // store the clipping frequency here
    uint64_t* clip1frequency = ckallocz(maxreadlength * sizeof(uint64_t));
    uint64_t* clip2frequency = ckallocz(maxreadlength * sizeof(uint64_t));

    int ret;
    uint64_t numread = 0;
    while((ret = bam_read1(fp, alignment)) >= 0){
        numread += 1;
        if((numread % 10000000) == 0){
            fprintf(stderr, "Processed %"PRIu64" reads\n", numread);
        }

        // if this is a zero length read, then I have nothing to do
        if(alignment->core.l_qseq == 0) continue;

        // do we only want the stats for some particular run?
        if(readstarts != NULL){
            if(0 != strncmp(bam1_qname(alignment), 
                            readstarts, 
                            strlen(readstarts))){
                continue;    
            }
        }

        // do we only want the stats for some particular read group?
        if(readgroups != NULL){
            uint8_t* rg = bam_aux_get(alignment, "RG");
            if(rg == NULL){
                fatal("Option -x specified & no read group for a read");
            }   
            forceassert(rg[0] == 'Z');
            if(0 != strcmp(bam_aux2Z(rg), readgroups)){
                continue;
            }
        }

        // if this is a secondary or supplementary alignment, then I do not want
        // to record this in the stats
        if(((alignment->core.flag & 0x100) == 0x100) || 
           ((alignment->core.flag & 0x800) == 0x800)){
            continue;
        }

        // if this is a QC failed read, just register that count and move on
        if((alignment->core.flag & 0x200) == 0x200){
            stats->r_qcfailed += 1;
            stats->b_qcfailed += alignment->core.l_qseq;
            continue;   
        }

        if(1 == debug_flag){
            fprintf(stderr, "Read name : %s\n", bam1_qname(alignment));
        }

        stats->r_qcpassed += 1;
        stats->b_qcpassed += alignment->core.l_qseq;

        // update the data structures to hold this read information
        if(alignment->core.l_qseq > maxreadlength){
            readlen1frequency = ckrealloc(readlen1frequency,
                                   alignment->core.l_qseq * sizeof(uint64_t));
            for(i = maxreadlength; i < alignment->core.l_qseq; i++){
                readlen1frequency[i] = 0;
            }
            readlen2frequency = ckrealloc(readlen2frequency,
                                   alignment->core.l_qseq * sizeof(uint64_t));
            for(i = maxreadlength; i < alignment->core.l_qseq; i++){
                readlen2frequency[i] = 0;
            }
            qual1values = ckrealloc(qual1values, 
                                    alignment->core.l_qseq * sizeof(uint64_t*));
            qual2values = ckrealloc(qual2values, 
                                    alignment->core.l_qseq * sizeof(uint64_t*));
            for(i = maxreadlength; i < alignment->core.l_qseq; i++){
                qual1values[i] = ckallocz(maxqualityvalue * sizeof(uint64_t));
                qual2values[i] = ckallocz(maxqualityvalue * sizeof(uint64_t));
            }
            nuc1counts = ckrealloc(nuc1counts,
                                   alignment->core.l_qseq * sizeof(uint64_t*));
            nuc2counts = ckrealloc(nuc2counts,
                                   alignment->core.l_qseq * sizeof(uint64_t*));
            for(i = maxreadlength; i < alignment->core.l_qseq; i++){
                nuc1counts[i] = ckallocz(5 * sizeof(uint64_t));
                nuc2counts[i] = ckallocz(5 * sizeof(uint64_t));
            }
            clip1frequency = ckrealloc(clip1frequency,
                                   alignment->core.l_qseq * sizeof(uint64_t));
            for(i = maxreadlength; i < alignment->core.l_qseq; i++){
                clip1frequency[i] = 0;
            }
            clip2frequency = ckrealloc(clip2frequency,
                                   alignment->core.l_qseq * sizeof(uint64_t));
            for(i = maxreadlength; i < alignment->core.l_qseq; i++){
                clip2frequency[i] = 0;
            }

            maxreadlength = alignment->core.l_qseq;
        }

        // is this a paired-end read or is it a single fragment?
        if((alignment->core.flag & 0x1) == 0){
            stats->r_singletons += 1;
            stats->b_singletons += alignment->core.l_qseq;
            
            // add the read length distribution
            forceassert(maxreadlength >= alignment->core.l_qseq);
            readlen1frequency[alignment->core.l_qseq - 1]++;

            // quality value distribution
            add_qualitycounts(alignment, qual1values);

            // nucleotide counts
            add_nucleotidecounts(alignment, nuc1counts);

            // clipping frequency
            add_clippingcounts(alignment, clip1frequency);
        }else{
            stats->r_paired += 1;
            stats->b_paired += alignment->core.l_qseq;

            // which fragment is it? read1 or read2
            if((alignment->core.flag & 0x40) == 0x40){
                assert((alignment->core.flag & 0x80) != 0x80);
                stats->read1s++;
        
                // add the read length distribution
                readlen1frequency[alignment->core.l_qseq - 1]++;

                // quality value distribution
                add_qualitycounts(alignment, qual1values);

                // nucleotide counts
                add_nucleotidecounts(alignment, nuc1counts);

                // clipping frequency
                add_clippingcounts(alignment, clip1frequency);
            }else if((alignment->core.flag & 0x80) == 0x80){
                assert((alignment->core.flag & 0x40) != 0x40);
                stats->read2s++;

                // add the read length distribution
                readlen2frequency[alignment->core.l_qseq - 1]++;

                // quality value distribution
                add_qualitycounts(alignment, qual2values);

                // nucleotide counts
                add_nucleotidecounts(alignment, nuc2counts);

                // clipping frequency
                add_clippingcounts(alignment, clip2frequency);
            }else{
                fatalf("A paired read is expected to have just 2 fragments\n");
            }
        }

        // nothing else to do if this did not align to the reference sequence
        if((alignment->core.flag & 0x4) == 0x4) continue;
        stats->r_aligned += 1;

        // how many bases in this read were used in this alignment?
        int nbases_used = lookatcigar(alignment);
        stats->b_aligned += nbases_used;

        // did it map to the reference? mapping refers to an unique alignment
        if(alignment->core.qual > 0 && alignment->core.qual != 255){
            stats->r_mapped += 1;
            stats->b_mapped += nbases_used;
        }

        // is this read a duplicate
        if((alignment->core.flag & 0x400) == 0x400){
            stats->r_duplicates += 1;
            stats->b_duplicates += nbases_used;
        }

        // is this part of a pair?
        if((alignment->core.flag & 0x1) == 0x1){
            // is this read properly paired?
            if((alignment->core.flag & 0x2) == 0x2){
                stats->properlypaired++;
                stats->bothaligned++;

                // if this is properly paired then lets store the value for the
                // insert length distribution
                if(alignment->core.isize > 0){
                    if(alignment->core.isize > maxinsertsize){
                        insertlenfrequency = ckrealloc(insertlenfrequency,
                                      alignment->core.isize * sizeof(uint64_t));
                        for(i = maxinsertsize; i < alignment->core.isize; i++){
                            insertlenfrequency[i] = 0;
                        }
                        maxinsertsize = alignment->core.isize;
                    }
                    insertlenfrequency[alignment->core.isize - 1] += 1;
                }
            }else{
                // both reads from this fragment align?
                if((alignment->core.flag & 0x8) == 0x0){
                    stats->bothaligned++;

                    if((allinsertlengths ==TRUE) && (alignment->core.isize >0)){
                        if(alignment->core.isize > maxinsertsize){
                            insertlenfrequency = ckrealloc(insertlenfrequency,
                                      alignment->core.isize * sizeof(uint64_t));
                            for(i = maxinsertsize;i < alignment->core.isize;i++){
                            insertlenfrequency[i] = 0;
                        }
                            maxinsertsize = alignment->core.isize;
                        }
                        insertlenfrequency[alignment->core.isize - 1] += 1;
                    }
                }
            }
        }

        // register this reads contribution to the coverage distribution; only
        // if this is not a duplicate read
        if(((alignment->core.flag & 0x400) != 0x400) && (reference != NULL)){
            int pos1 = alignment->core.pos;
            int pos2 = bam_calend(&alignment->core, bam1_cigar(alignment)); 
            
            chrcoverage* cov = must_find_hashtable(reference,
                               hin->target_name[alignment->core.tid],
                               strlen(hin->target_name[alignment->core.tid]));

            for(i = pos1; i < pos2; i++){
                if(i < (int)cov->length){
                    cov->coverage[i] += 1;            
                    if(cov->coverage[i] == 251) cov->coverage[i] = 250;
                }
            }
        }
    }

    bam_close(fp);
    //timestamp("Read the BAM file\n");

    // print the stats for the BAM file
    if(statsname == NULL){
        print_stats(stats, stdout);
    }else{
        FILE* fp = ckopen(statsname, "w");
        print_stats(stats, fp);
        fclose(fp);
    }

    // print the read length distribution 
    print_rlength(rlensname, 
                  readlen1frequency, 
                  readlen2frequency, 
                  maxreadlength);

    // print the quality value distribution
    print_readqualdist(qualsname, qual1values, qual2values, maxreadlength);

    // print the nucleotide composition
    print_nuccomp(nucsname, nuc1counts, nuc2counts, maxreadlength);

    // print the insert length distribution
    print_insertlendist(insname, insertlenfrequency, maxinsertsize);

    // print the clipping frequency
    print_clippingcounts(clipname, clip1frequency, clip2frequency,maxreadlength);

    if(reference != NULL){
        // print the coverage distribution over the reference sequence
        print_coverage_distribution(covname, reference);

        // print the binned GC coverage in non-overlapping windows
        print_gccoverage(gccovname, reference, windowsize);
    }
}

int main(int argc, char** argv)
{
    argv0 = "bamstats";
    int c;

    char* rlensname = NULL;
    char* qualsname = NULL;
    char* nucsname  = NULL;
    char* covname   = NULL;
    char* insname   = NULL;
    char* gccovname = NULL;
    char* statsname = NULL;
    char* clipname  = NULL;
    int windowsize  = 1000000;
    char* readstarts= NULL;
    char* readgroups = NULL;
    nowarnings = FALSE;

    // by default I only include the properly paired reads in the insert length
    // distribution. But in some cases (like damaged ancient reads), very few
    // reads get marked as "properly paired" since BWA does not find enough
    // pairs to calculate the mean and sd values that should be called as
    // properly paired. In those cases I would want to include all pairs in this
    // calculation. Of course that means all pairs where both ends mapped.
    bool allinsertlengths = FALSE;

    while (1){
        static struct option long_options[] = {
            {"debug"   , no_argument , 0, 'd'},
            {"rlens"   , required_argument, 0, 'l'},
            {"rqual"   , required_argument, 0, 'q'},
            {"rnucs"   , required_argument, 0, 'n'},
            {"rcovs"   , required_argument, 0, 'c'},
            {"rinst"   , required_argument, 0, 'i'},
            {"gccov"   , required_argument, 0, 'g'},  
            {"stats"   , required_argument, 0, 's'},
            {"wsize"   , required_argument, 0, 'w'},
            {"run"     , required_argument, 0, 'r'},
            {"rallinst", no_argument,       0, 'a'},
            {"nowarn"  , no_argument,       0, 'f'},
            {"rg"      , required_argument, 0, 'x'},
            {"clip"    , required_argument, 0, 'y'},
            {0, 0, 0, 0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "dl:q:n:c:i:g:s:w:r:afx:y:",
                        long_options, &option_index);

        if (c == -1) break;

        switch (c){
            case 0:
                break;
            case 'd':
                debug_flag = TRUE;
                break;
            case 'l':
                rlensname  = optarg;
                break;
            case 'q':
                qualsname  = optarg;
                break;
            case 'n':
                nucsname   = optarg;
                break;
            case 'c':
                covname    = optarg;
                break;
            case 'i':
                insname    = optarg;
                break;
            case 'g':
                gccovname  = optarg;
                break;
            case 's':
                statsname  = optarg;
                break;
            case 'w':
                windowsize = atoi(optarg);
                break;
            case 'r':
                readstarts = optarg;
                break;
            case 'y':
                clipname = optarg;
                break;
            case 'a':
                allinsertlengths = TRUE;
                break;
            case 'f':
                nowarnings = TRUE;
                break;
            case 'x':
                readgroups = optarg;
                break;
            case '?':
                break;
            default:
                abort();
        }
    }

    /* start clock book-keeping */
    t0 = time(0);

    if((argc - optind) != 2){
        fprintf(stderr, "Usage\n");
        fprintf(stderr, "\tbamstats [options] reference.fa alignments.bam\n");
        return EXIT_FAILURE;
    }

    // user should have only specified one of readgroups or readstarts
    if((readgroups != NULL ) && (readstarts != NULL)){
        fprintf(stderr, "This will count reads which have the given startname AND the given readgroup\n");
    }

    calcbamstats(argv[optind], 
                 argv[optind+1], 
                 rlensname,
                 qualsname, 
                 nucsname, 
                 covname, 
                 insname,
                 gccovname,
                 clipname,
                 windowsize,
                 statsname,
                 readstarts,
                 readgroups,
                 allinsertlengths);

    /* print the relevant stats used by the program */
    print_usage();
   
    return EXIT_SUCCESS;
}
