// Differences from samtools flagstat
// a) Supplementary and secondary alignments are not ignored by samtools
//    flagstat when calculating QC-passed reads.
// b) samtools mapped is the same as the aligned as output by this tool. mapped
//    here refers to alignments with non-zero mapping quality.
// Other notes:
// a) The nucleotide and quality distribution is only made for reads that pass
//    qc.
// b) The quality values vary from 0 to 60.

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include <locale.h>

#include "sam.h"

#include "errors.h"
#include "time.h"
#include "memalloc.h"
#include "files.h"
#include "sequences.h"
#include "hashtable.h"

#define MAXQ 60

// print debug information?
int debug_flag = FALSE;

typedef int64_t i64;
typedef int32_t i32;
typedef uint64_t u64;
typedef uint32_t u32;

typedef struct bamstats_t{
    // generation statistics
    u64 r_qcpassed;
    u64 r_qcfailed;
    u64 r_paired;
    u64 r_singletons;
    u64 b_qcpassed;
    u64 b_qcfailed;
    u64 b_paired;
    u64 b_singletons;
    
    // alignment statistics
    u64 r_duplicates;
    u64 r_aligned;
    u64 r_mapped;
    u64 b_duplicates;
    u64 b_aligned;
    u64 b_mapped;

    u64 read1s;
    u64 read2s;
    u64 properlypaired;
    u64 bothaligned;
}bamstats;

typedef struct dist_st{
    u64* r1length;
    u64* r2length;  
    u64* insrtlength;
    u64** n1counts;
    u64** n2counts;
    u64** q1counts;
    u64** q2counts;
    u64* c1counts;
    u64* c2counts;
    u64* m1counts;
    u64* m2counts;
    u64* i1counts;  
    u64* i2counts;
    i32 max_insert_size;
    i32 max_read_length;
}dist;

typedef struct chrcoverage_t
{
    uint64_t length;
    uchar* coverage;
    uchar* sequence;
} chrcoverage;

typedef struct wflags_st{
    bool q_past_threshold;
    bool noseq;
    bool noqual;
}wflags;
wflags* warnings = NULL;

// if locale is not supported, then I could use this function. 
void printfcomma (int n) {
    int n2 = 0;
    int scale = 1;
    if (n < 0) {
        printf ("-");
        n = -n;
    }
    while (n >= 1000) {
        n2 = n2 + scale * (n % 1000);
        n /= 1000;
        scale *= 1000;
    }
    printf ("%d", n);
    while (scale != 1) {
        scale /= 1000;
        n = n2 / scale;
        n2 = n2  % scale;
        printf (",%03d", n);
    }
}

// print the statistics for this BAM file
static void print_stats(const bamstats* stats, FILE* fp)
{
    setlocale(LC_NUMERIC, "");
    fprintf(fp, "Reads Summary\n");
    fprintf(fp,"-------------\n");
    fprintf(fp,"QC-passed reads : %'"PRIu64"\n", stats->r_qcpassed);
    fprintf(fp,"QC-failed reads : %'"PRIu64"\n", stats->r_qcfailed);
    fprintf(fp,"Number of single-fragments : %'"PRIu64"\n", stats->r_singletons);
    fprintf(fp,"Number of reads paired in sequencing : %'"PRIu64"\n", stats->r_paired);
    fprintf(fp,"\n");
    fprintf(fp,"Duplicate reads : %'"PRIu64" (%2.2f%%)\n", 
            stats->r_duplicates, 
            stats->r_duplicates * 100.0 / stats->r_qcpassed);
    fprintf(fp,"Number of reads aligned : %'"PRIu64" (%2.2f%%)\n", 
            stats->r_aligned,
            stats->r_aligned * 100.0 / stats->r_qcpassed);
    fprintf(fp,"Number of reads mapped  : %'"PRIu64" (%2.2f%%)\n", 
            stats->r_mapped,
            stats->r_mapped * 100.0 / stats->r_qcpassed);
    fprintf(fp,"\n");
    fprintf(fp,"Number of read1's : %'"PRIu64"\n", stats->read1s);
    fprintf(fp,"Number of read2's : %'"PRIu64"\n", stats->read2s);
    fprintf(fp,"Number of reads properly paired : %'"PRIu64" (%2.2f%%)\n", 
            stats->properlypaired,
            stats->properlypaired * 100.0 / stats->r_paired);
    fprintf(fp,"Number of reads with itself and mate aligned : %'"PRIu64" (%2.2f%%)\n", 
            stats->bothaligned,
            stats->bothaligned * 100.0 / stats->r_paired);
    fprintf(fp,"\n");
    fprintf(fp,"Bases Summary\n");
    fprintf(fp,"-------------\n");
    fprintf(fp,"QC-passed bases : %'"PRIu64"\n", stats->b_qcpassed);
    fprintf(fp,"QC-failed bases : %'"PRIu64"\n", stats->b_qcfailed);
    fprintf(fp,"Number of bases in single-fragments : %'"PRIu64"\n", stats->b_singletons);
    fprintf(fp,"Number of bases paired in sequencing : %'"PRIu64"\n", stats->b_paired);
    fprintf(fp,"\n");
    fprintf(fp,"Duplicate bases : %'"PRIu64" (%2.2f%%)\n", 
            stats->b_duplicates, 
            stats->b_duplicates * 100.0 / stats->b_qcpassed);
    fprintf(fp,"Number of bases aligned : %'"PRIu64" (%2.2f%%)\n", 
            stats->b_aligned,
            stats->b_aligned * 100.0 / stats->b_qcpassed);
    fprintf(fp,"Number of bases mapped  : %'"PRIu64" (%2.2f%%)\n", 
            stats->b_mapped,
            stats->b_mapped * 100.0 / stats->b_qcpassed);
    fprintf(fp,"\n");
}

// return the number of bases in the read that are used in the alignment
static int lookatcigar(const bam1_t* const alignment)
{
    int numused = 0;

    u32* cigar = bam_get_cigar(alignment);

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

static void add_quality_counts(const bam1_t* const alignment,
                               dist* const distributions)
{
    uchar* qual = bam_get_qual(alignment);

    u64** qualcounts = NULL;
    if ((alignment->core.flag & 0x80) == 0x80) {
        qualcounts = distributions->q2counts;
    } else {
        qualcounts = distributions->q1counts;
    }

    if ((alignment->core.flag & 0x10) == 0x10) {
        // the query is reverse complemented
        int i, j;
        for(i = alignment->core.l_qseq - 1, j = 0; i >= 0; i--, j++){
            if((warnings->q_past_threshold == FALSE) && (qual[i] >= MAXQ)) {
                warnf("Quality values in this BAM exceed threshold of 60\n");
                warnings->q_past_threshold = TRUE;
            }
            qualcounts[j][qual[i] >= MAXQ ? MAXQ-1 : qual[i]] += 1;
        }
    } else {
        int i;
        for(i = 0; i < alignment->core.l_qseq; i++){
            if((warnings->q_past_threshold == FALSE) && (qual[i] >= MAXQ)) {
                warnf("Quality values in this BAM exceed threshold of 60\n");
                warnings->q_past_threshold = TRUE;
            }
            qualcounts[i][qual[i] >= MAXQ ? MAXQ-1 : qual[i]] += 1;
        }
    }
}

static void add_clipping_locations(const bam1_t* const alignment,
                                   dist* const distributions)
{
    u32* cigar = NULL;
    u32 op, oplen;

    cigar = bam_get_cigar(alignment);
    if ((alignment->core.flag & 0x10) == 0x10) {
        op = cigar[alignment->core.n_cigar-1] & 0xf;
        oplen = (cigar[alignment->core.n_cigar-1] & 0xfffffff0) >> 4;
    } else {
        op = cigar[0] & 0xf;
        oplen = (cigar[0] & 0xfffffff0) >> 4;
    }

    if (op == BAM_CSOFT_CLIP) {
        if ((alignment->core.flag & 0x80) == 0x80) {
            distributions->c2counts[oplen] += 1;
        } else {
            distributions->c1counts[oplen] += 1;
        }  
    }
}

static inline int encodedbase_to_ascii(const int base) 
{
    if (base == 1) {
        return 'A';
    } else if (base == 2) {
        return 'C';
    } else if (base == 4) {
        return 'G';
    } else if (base == 8) {
        return 'T';
    } else if (base == 15) {
        return 'N';
    } else {
        fatalf("Unknown encoded base");
    }
    return 0;
}

static void add_indel_counts(const bam1_t* const alignment,
                            dist* const distributions)
{
    uchar* qseq = bam_get_seq(alignment);
    if ((*qseq == '*') && (warnings->noseq == FALSE)){
        warnf("The read sequences is not stored in this BAM for this read %s\n",
               bam_get_qname(alignment));
        warnings->noseq = TRUE;
    }

    u64* errfrequency = NULL;
    if ((alignment->core.flag & 0x80) == 0x80) {
        errfrequency = distributions->i2counts;
    } else {    
        errfrequency = distributions->i1counts;
    }           
    u32* cigar = bam_get_cigar(alignment);

    bool is_reverse = (alignment->core.flag & 0x10) == 0x10;

    uint i,oplen;
    uint qindx = 0;
    uint rindx = alignment->core.pos;
    for(i = 0; i < alignment->core.n_cigar; i++){
        oplen = ((cigar[i] & 0xfffffff0) >> 4);
        switch(cigar[i] & 0xf){
            case BAM_CMATCH :
                qindx += oplen;
                rindx += oplen;
                break;
            case BAM_CINS:
                if (is_reverse == FALSE) {
                    errfrequency[qindx] += 1;
                } else {
                    errfrequency[alignment->core.l_qseq - qindx-1] += 1;
                }
                qindx += oplen; 
                break;
            case BAM_CDEL:
                if (is_reverse == FALSE) {
                    errfrequency[qindx] += 1;
                } else {
                    errfrequency[alignment->core.l_qseq - qindx-1] += 1;
                }
                rindx += oplen; 
                break;
            case BAM_CREF_SKIP:
                fatalf("Have not handled REF_SKIP in CIGAR");
            case BAM_CSOFT_CLIP:
                qindx += oplen;
                break;
            case BAM_CHARD_CLIP:
                qindx += oplen;
                break;
            case BAM_CPAD:
                break;
            default :
                fatalf("Unhandled CIGAR operation\n");
        }
    }

}

static void add_mismatch_counts(const bam1_t* const alignment,
                                hashtable* const reference,
                                const bam_hdr_t* const hin,
                                dist* const distributions)
{
    chrcoverage* cov = must_find_hashtable(reference,
                               hin->target_name[alignment->core.tid],
                               strlen(hin->target_name[alignment->core.tid]));
    uchar* qseq = bam_get_seq(alignment);
    if ((*qseq == '*') && (warnings->noseq == FALSE)){
        warnf("The read sequences is not stored in this BAM for this read %s\n",
               bam_get_qname(alignment));
        warnings->noseq = TRUE;
    }

    u64* errfrequency = NULL;
    if ((alignment->core.flag & 0x80) == 0x80) {
        errfrequency = distributions->m2counts;
    } else {    
        errfrequency = distributions->m1counts;
    }           
    u32* cigar = bam_get_cigar(alignment);

    bool is_reverse = (alignment->core.flag & 0x10) == 0x10;

    uint i, j, oplen;
    uint qindx = 0;
    uint rindx = alignment->core.pos;
    for(i = 0; i < alignment->core.n_cigar; i++){
        oplen = ((cigar[i] & 0xfffffff0) >> 4);
        switch(cigar[i] & 0xf){
            case BAM_CMATCH :
                for (j = 0; j < oplen; j++) {
                    int q = encodedbase_to_ascii(bam_seqi(qseq,qindx));
                    int r = toupper(cov->sequence[rindx]);
                    if (q != r){
                        if (is_reverse == FALSE) {
                            errfrequency[qindx] += 1;
                        } else {
                            errfrequency[alignment->core.l_qseq - qindx-1] += 1;
                        }
                    }
                    qindx += 1;
                    rindx += 1;
                }
                break;
            case BAM_CINS:
                qindx += oplen; 
                break;
            case BAM_CDEL:
                rindx += oplen; 
                break;
            case BAM_CREF_SKIP:
                fatalf("Have not handled REF_SKIP in CIGAR");
            case BAM_CSOFT_CLIP:
                qindx += oplen;
                break;
            case BAM_CHARD_CLIP:
                qindx += oplen;
                break;
            case BAM_CPAD:
                break;
            default :
                fatalf("Unhandled CIGAR operation\n");
        }
    }
}

static void add_nucleotide_counts(const bam1_t* const alignment, 
                                  dist* const distributions) 
{
    uchar* seq = bam_get_seq(alignment);
    if((*seq == '*') && (warnings->noseq == FALSE)){
        warnf("The read sequences is not stored in this BAM for this read %s\n",
                bam_get_qname(alignment));
        warnings->noseq = TRUE;
    }
    if(*seq == '*') return;

    u64** nuccounts = NULL;
    if ((alignment->core.flag & 0x80) == 0x80) {
        nuccounts = distributions->n2counts;
    } else {
        nuccounts = distributions->n1counts;
    }

    if ((alignment->core.flag & 0x10) == 0x10) {
        // the query is reverse complemented
        int i, j;
        for(i = alignment->core.l_qseq - 1, j = 0; i >= 0; i--, j++){
            switch(bam_seqi(seq, i)){
                case 1 : nuccounts[j][0] += 1; break;
                case 2 : nuccounts[j][1] += 1; break;
                case 4 : nuccounts[j][2] += 1; break;
                case 8 : nuccounts[j][3] += 1; break;
                case 15: nuccounts[j][4] += 1; break;
                default: 
                    fatalf("Incorrect base call.\n");
            }   
        }
    } else {
        int i;
        for(i = 0; i < alignment->core.l_qseq; i++){
            switch(bam_seqi(seq, i)){
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
}

static void print_read_lengths(const dist* const distributions,
                               const char* const rlensname)
{
    if (rlensname == NULL) return;

    FILE* fp = ckopen(rlensname, "w");
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        fprintf(fp, "1\t%u\t%"PRIu64"\n", i+1, distributions->r1length[i]);
    }
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        fprintf(fp, "2\t%u\t%"PRIu64"\n", i+1, distributions->r2length[i]);
    }
    fclose(fp);
}

static void print_insert_dist(const dist* const distributions,
                              const char* const insrtname) 
{
    if (insrtname == NULL) return;

    FILE* fp = ckopen(insrtname, "w");
    for (i32 i = 0; i < distributions->max_insert_size; i++) {
        fprintf(fp, "%u\t%"PRIu64"\n", i+1, distributions->insrtlength[i]);
    }
    fclose(fp);
}

static void print_nuc_distribution(const dist* const distributions, 
                                   const char* const rnucsname)
{
    if (rnucsname == NULL) return;

    FILE* fp = ckopen(rnucsname, "w");
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        fprintf(fp, 
     "1\t%u\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",
                i+1, 
                distributions->n1counts[i][0],
                distributions->n1counts[i][1],
                distributions->n1counts[i][2],
                distributions->n1counts[i][3],
                distributions->n1counts[i][4],
                distributions->n1counts[i][0] +
                distributions->n1counts[i][1] +
                distributions->n1counts[i][2] +
                distributions->n1counts[i][3] +
                distributions->n1counts[i][4]);
    }
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        fprintf(fp, 
     "2\t%u\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n",
                i+1, 
                distributions->n2counts[i][0],
                distributions->n2counts[i][1],
                distributions->n2counts[i][2],
                distributions->n2counts[i][3],
                distributions->n2counts[i][4],
                distributions->n2counts[i][0] +
                distributions->n2counts[i][1] +
                distributions->n2counts[i][2] +
                distributions->n2counts[i][3] +
                distributions->n2counts[i][4]);
    }
    fclose(fp);
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
    for(uint i = 0; i < MAXQ; i++){
        numentries += frequencytable[i];
    }
    if(numentries == 0) {
        *pwmin = 0;
        *pbmin = 0;
        *pmedian = 0;
        *pbmax = 0;
        *pwmax = 0;
        return;
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
    for(i = 0; i < (int)MAXQ; i++){
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

static void print_qual_distribution(const dist* const distributions,
                                    const char* const rqualname)
{
    if (rqualname == NULL) return;
    
    FILE* fp = ckopen(rqualname, "w");
    int whisker_min, box_min, median, box_high, whisker_high;
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        summary(distributions->q1counts[i],
                &whisker_min, &box_min, &median, &box_high, &whisker_high);
        fprintf(fp, "1\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    i+1, whisker_min, box_min, median, box_high, whisker_high);
        }
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        summary(distributions->q2counts[i],
                &whisker_min, &box_min, &median, &box_high, &whisker_high);
        fprintf(fp, "2\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    i+1, whisker_min, box_min, median, box_high, whisker_high);
        }
    fclose(fp);
}

// print the coverage distribution of the reference sequence
static void print_coverage_distribution(const hashtable* const reference,
                                        const char* const covname)
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
        fprintf(fp, "%d\t%u\n", i, coverage[i]);
    }

    fclose(fp);
    ckfree(coverage);
}

static void print_fclip_distribution(const dist* const distributions,
                                     const char* const fclipname)
{
    if (fclipname == NULL) return;

    FILE* fp = ckopen(fclipname, "w");
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        u64 total = distributions->n1counts[i][0] + distributions->n1counts[i][1]                  + distributions->n1counts[i][2] + distributions->n1counts[i][3]
                  + distributions->n1counts[i][4];
        fprintf(fp, "1\t%u\t%"PRIu64"\t%"PRIu64"\n", 
                i+1, distributions->c1counts[i],total);
    }
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        u64 total = distributions->n2counts[i][0] + distributions->n2counts[i][1]                  + distributions->n2counts[i][2] + distributions->n2counts[i][3]
                  + distributions->n2counts[i][4];
        fprintf(fp, "2\t%u\t%"PRIu64"\t%"PRIu64"\n", 
                    i+1, distributions->c2counts[i], total);
    }
    fclose(fp);
}

static void print_mm_distribution(const dist* const distributions, 
                                  const char* const msmtcname)
{
    if (msmtcname == NULL) return;

    FILE* fp = ckopen(msmtcname, "w");
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        u64 total = distributions->n1counts[i][0] + distributions->n1counts[i][1]                  + distributions->n1counts[i][2] + distributions->n1counts[i][3]
                  + distributions->n1counts[i][4];
        fprintf(fp, "1\t%u\t%"PRIu64"\t%"PRIu64"\n", 
                    i+1, distributions->m1counts[i], total);
    }
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        u64 total = distributions->n2counts[i][0] + distributions->n2counts[i][1]                  + distributions->n2counts[i][2] + distributions->n2counts[i][3]
                  + distributions->n2counts[i][4];
        fprintf(fp, "2\t%u\t%"PRIu64"\t%"PRIu64"\n", 
                    i+1, distributions->m2counts[i], total);
    }
    fclose(fp);
}

static void print_indel_distribution(const dist* const distributions,
                                  const char* const indelname)
{
    if (indelname == NULL) return;

    FILE* fp = ckopen(indelname, "w");
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        u64 total = distributions->n1counts[i][0] + distributions->n1counts[i][1]                  + distributions->n1counts[i][2] + distributions->n1counts[i][3]
                  + distributions->n1counts[i][4];
        fprintf(fp, "1\t%u\t%"PRIu64"\t%"PRIu64"\n", 
                    i+1, distributions->i1counts[i], total);
    }
    for (i32 i = 0; i < distributions->max_read_length; i++) {
        u64 total = distributions->n2counts[i][0] + distributions->n2counts[i][1]                  + distributions->n2counts[i][2] + distributions->n2counts[i][3]
                  + distributions->n2counts[i][4];
        fprintf(fp, "2\t%u\t%"PRIu64"\t%"PRIu64"\n", 
                    i+1, distributions->i2counts[i], total);
    }
    fclose(fp);
}

static void print_gc_coverage(const hashtable* const reference,
                              const char* const gccovname,
                              const u64 windowsize)
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

static hashtable* read_reference(const char* const refname)
{
    hashtable* reference = new_hashtable(12);

    sequence* sp = read_fasta_sequence(refname);
    
    while(sp != NULL){
        // allocate a coverage array for the sequence
        chrcoverage* cov = ckallocz(sizeof(chrcoverage));
        cov->length   = strlen((char*)sp->sequence);
        cov->coverage = ckallocz(strlen((char*)sp->sequence));
        cov->sequence = ckallocz(strlen((char*)sp->sequence)+1);
        cov->sequence = (uchar*)strcpy((char*)cov->sequence,(char*)sp->sequence);

        // if the name of the sequence has more than one tokens, just use the
        // first token in the name
        int i = 0;
        while((sp->header[i] != '\n') && 
              (sp->header[i] != 0)    && 
              (sp->header[i] != '\t') && 
              (sp->header[i] != 32)) i++;
        sp->header[i] = 0;

        add_hashtable(reference,(char*)sp->header,strlen((char*)sp->header),cov);
        sp = get_next_sequence(sp);
    } 

    return reference;
}

static int calcbamstats(const char* const refname,
                        const char* const bamname,
                        const char* const statsname,
                        const char* const rlensname,
                        const char* const insrtname,
                        const char* const rnucsname,
                        const char* const rqualname,
                        const char* const rcovsname,
                        const char* const gccovname,
                        const u64 windowsize,
                        const char* const fclipname,
                        const char* const msmtcname,
                        const char* const indelname)
{
    // lets read the reference sequence. The reference sequence will be a
    // hashtable where the key is the name of the chromosome and the value will
    // be a char array of size equal to  the length of the chromosome sequence.
    hashtable* reference = NULL;
    if((rcovsname != NULL) || (gccovname != NULL) || (msmtcname != NULL)){
        reference = read_reference(refname);
    }

    int status;
    samFile *in = sam_open(bamname, "r");
    bam_hdr_t *hdr = NULL;
    if (in == NULL) {
        status = EXIT_FAILURE;
        return 0;
    }
    hdr = sam_hdr_read(in);
    if (hdr == NULL) {
        status = EXIT_FAILURE;
        goto clean;
    }

    int ret;
    bam1_t *b = bam_init1();
    u64 numread = 0; // number of reads analyzed
    bamstats* stats = ckallocz(sizeof(bamstats)); // stats
    dist* distributions = ckallocz(sizeof(dist)); // distributions

    while ((ret = sam_read1(in, hdr, b)) >= 0) {
        numread += 1;
        if ((numread % 1000000) == 0) {
            fprintf(stderr, "Processed %"PRIu64" reads\n", numread);
        }
        if (1 == debug_flag) {
            fprintf(stderr, "Read name : %s\n", bam_get_qname(b));
        }

        // ignore if this is a zero length read (have seen it in some cases
        // where the reads were clipped by another tool. also ignore all
        // secondary or supplementary alignments for now. 
        // * samtools flagstat does not discard these (secondary/supp).
        if (b->core.l_qseq == 0) continue;
        if (((b->core.flag & 0x100) == 0x100) || 
            ((b->core.flag & 0x800) == 0x800)){
            continue;
        }
        
        // if this is a QC failed read, just register that count and move on
        if ((b->core.flag & 0x200) == 0x200) {
            stats->r_qcfailed += 1;
            stats->b_qcfailed += b->core.l_qseq;
            continue;   
        }
        stats->r_qcpassed += 1;
        stats->b_qcpassed += b->core.l_qseq;

        if (b->core.l_qseq > distributions->max_read_length) {
            distributions->r1length = ckreallocz(distributions->r1length,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64)); 
            distributions->r2length = ckreallocz(distributions->r2length,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));

            distributions->n1counts = ckreallocz(distributions->n1counts,
                                  distributions->max_read_length * sizeof(u64*),
                                                  b->core.l_qseq * sizeof(u64*));
            distributions->n2counts = ckreallocz(distributions->n2counts,
                                  distributions->max_read_length * sizeof(u64*),
                                                  b->core.l_qseq * sizeof(u64*));
            for (i32 i = distributions->max_read_length;
                     i < b->core.l_qseq; i++) {
                distributions->n1counts[i] = ckallocz(5 * sizeof(u64));
                distributions->n2counts[i] = ckallocz(5 * sizeof(u64));
            }
         
            distributions->q1counts = ckreallocz(distributions->q1counts,
                                  distributions->max_read_length * sizeof(u64*),
                                                  b->core.l_qseq * sizeof(u64*));
            distributions->q2counts = ckreallocz(distributions->q2counts,
                                  distributions->max_read_length * sizeof(u64*),
                                                  b->core.l_qseq * sizeof(u64*));
            for (i32 i = distributions->max_read_length;
                     i < b->core.l_qseq; i++) {
                distributions->q1counts[i] = ckallocz(MAXQ * sizeof(u64));
                distributions->q2counts[i] = ckallocz(MAXQ * sizeof(u64));
            }

            distributions->c1counts = ckreallocz(distributions->c1counts,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));
            distributions->c2counts = ckreallocz(distributions->c2counts,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));
   
            distributions->m1counts = ckreallocz(distributions->m1counts,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));
            distributions->m2counts = ckreallocz(distributions->m2counts,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));
 
            distributions->i1counts = ckreallocz(distributions->i1counts,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));
            distributions->i2counts = ckreallocz(distributions->i2counts,
                                 distributions->max_read_length * sizeof(u64),
                                                 b->core.l_qseq * sizeof(u64));
       
            distributions->max_read_length = b->core.l_qseq;
        }

        // add the nucleotide counts from this read
        add_nucleotide_counts(b, distributions);

        // add the quality distribution information
        add_quality_counts(b, distributions);

        // is this a paired-end read or is it a single fragment?
        if ((b->core.flag & 0x1) == 0) {
            stats->r_singletons += 1;
            stats->b_singletons += b->core.l_qseq;
            stats->read1s++;
            distributions->r1length[b->core.l_qseq-1] += 1;
        } else {
            stats->r_paired += 1;
            stats->b_paired += b->core.l_qseq;
            if ((b->core.flag & 0x40) == 0x40) {
                stats->read1s++;
                distributions->r1length[b->core.l_qseq-1] += 1;
            } else if ((b->core.flag & 0x80) == 0x80) {
                stats->read2s++;
                distributions->r2length[b->core.l_qseq-1] += 1;
            }

            // is this properly paired?
            if((b->core.flag & 0x2) == 0x2){
                stats->properlypaired++;
                stats->bothaligned++;

                if (b->core.isize >= distributions->max_insert_size) {
                    distributions->insrtlength =
                                    ckreallocz(distributions->insrtlength, 
                                  distributions->max_insert_size * sizeof(u64), 
                                                   b->core.isize * sizeof(u64));
                    distributions->max_insert_size = b->core.isize;
                }
                if (b->core.isize > 0) 
                    distributions->insrtlength[b->core.isize-1] += 1;
            }else{
                // both reads from this fragment align?
                if((b->core.flag & 0x8) == 0x0){
                    stats->bothaligned++;
                }
            }
        }

        // nothing else to do if this did not align to the reference sequence
        if((b->core.flag & 0x4) == 0x4) continue;
        stats->r_aligned += 1;

        // how many bases in this read were used in this alignment?
        int nbases_used = lookatcigar(b);
        stats->b_aligned += nbases_used;

        // did it map to the reference? mapping refers to an unique alignment
        if(b->core.qual > 0 && b->core.qual != 255){
            stats->r_mapped += 1;
            stats->b_mapped += nbases_used;
        }

        // is this read a duplicate
        if((b->core.flag & 0x400) == 0x400){
            stats->r_duplicates += 1;
            stats->b_duplicates += nbases_used;
            continue;
        }

        // add the 5' clipping locations on this read which is aligned
        add_clipping_locations(b, distributions);

        if (reference != NULL) {
            // register this reads contribution towards the coverage
            i32 pos1 = b->core.pos;
            i32 pos2 = bam_endpos(b); 
            chrcoverage* cov = must_find_hashtable(reference,
                               hdr->target_name[b->core.tid],
                               strlen(hdr->target_name[b->core.tid]));
            for(i32 i = pos1; i < pos2; i++){
                if(i < (int)cov->length){
                    cov->coverage[i] += 1;            
                    if(cov->coverage[i] == 251) cov->coverage[i] = 250;
                }
            }
            
            // add the mismatches
            add_mismatch_counts(b, reference, hdr, distributions);
        }
        // add any indels is present
        add_indel_counts(b, distributions);

    }
    bam_destroy1(b);
    if (ret != -1) // eof
        status = EXIT_FAILURE;

    // print the stats
    if(statsname == NULL){ 
        print_stats(stats, stdout);
    }else{
        FILE* fp = ckopen(statsname, "w");
        print_stats(stats, fp);
        fclose(fp);
    }
    
    // print the read length distribution
    print_read_lengths(distributions, rlensname);

    // print the insert length distribution
    print_insert_dist(distributions, insrtname);

    // print the nucleotide distribution
    print_nuc_distribution(distributions, rnucsname);

    // print the quality distribution
    print_qual_distribution(distributions, rqualname);

    // print the clipping information
    print_fclip_distribution(distributions, fclipname);

    if (reference != NULL) {
        // print the coverage distribution
        print_coverage_distribution(reference, rcovsname);

        // print the GC coverage information
        print_gc_coverage(reference, gccovname, windowsize);

        // print the mismatch location distribution
        print_mm_distribution(distributions, msmtcname);
    }
    // print the indel count distribution
    print_indel_distribution(distributions, indelname);


 clean:
    if (hdr != NULL) bam_hdr_destroy(hdr);
    if (hts_close(in) != 0)
        status = EXIT_FAILURE;
    return 1;

}

int main(int argc, char** argv)
{
    argv0 = "bamstats";
    int c;

    char* statsname = NULL;
    char* rlensname = NULL;
    char* insrtname = NULL;
    char* rnucsname = NULL;
    char* rqualname = NULL;
    char* rcovsname = NULL;
    char* gccovname = NULL;
    char* fclipname = NULL;
    char* msmtcname = NULL;
    char* indelname = NULL;
    u64 windowsize  = 1000000;

    while (1){
        static struct option long_options[] = {
            {"debug"   , no_argument , 0, 'd'},
            {"stats"   , required_argument, 0, 's'},
            {"rlens"   , required_argument, 0, 'r'},
            {"rinst"   , required_argument, 0, 'i'},
            {"rnucs"   , required_argument, 0, 'n'},
            {"rqual"   , required_argument, 0, 'q'},
            {"rcovs"   , required_argument, 0, 'c'},
            {"gccov"   , required_argument, 0, 'g'},
            {"wsize"   , required_argument, 0, 'w'},
            {"flcip"   , required_argument, 0, 'f'},
            {"mmtch"   , required_argument, 0, 'm'},
            {"indel"   , required_argument, 0, 'x'},
            {0, 0, 0, 0}
        };
    
        int option_index = 0;
        c = getopt_long(argc, argv, "ds:r:i:n:q:c:g:w:f:m:x:",
                        long_options, &option_index);

        if (c == -1) break;

        switch (c){
            case 0:
                break;
            case 'd':
                debug_flag = TRUE;
                break;
            case 's':
                statsname  = optarg;
                break;
            case 'r':
                rlensname = optarg;
                break;
            case 'i':
                insrtname = optarg;
                break;
            case 'n':
                rnucsname = optarg;
                break;
            case 'q':
                rqualname = optarg;
                break;
            case 'c':
                rcovsname = optarg;
                break;
            case 'g':
                gccovname = optarg;
                break;
            case 'm':
                msmtcname = optarg;
                break;
            case 'x':
                indelname = optarg;
                break;
            case 'w':
                if (sscanf(optarg, "%"PRIu64, &windowsize) != 1) {
                    fatalf("could not parse the windowsize: %s", optarg);
                }
                break;
            case 'f':
                fclipname = optarg;
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

    /* warnings on */
    warnings = ckallocz(sizeof(warnf));

    calcbamstats(argv[optind], 
                 argv[optind+1], 
                 statsname,
                 rlensname,
                 insrtname,
                 rnucsname,
                 rqualname,
                 rcovsname,
                 gccovname,
                 windowsize,
                 fclipname,
                 msmtcname,
                 indelname);

    /* print the relevant stats used by the program */
    print_usage();
   
    ckfree(warnings);
    return EXIT_SUCCESS;
}
