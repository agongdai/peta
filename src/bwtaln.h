#ifndef BWTALN_H
#define BWTALN_H

#include <stdint.h>

#define BWA_TYPE_NO_MATCH 0
#define BWA_TYPE_UNIQUE 1
#define BWA_TYPE_REPEAT 2
#define BWA_TYPE_MATESW 3

#define SAM_FPD   1 // paired
#define SAM_FPP   2 // properly paired
#define SAM_FSU   4 // self-unmapped
#define SAM_FMU   8 // mate-unmapped
#define SAM_FSR  16 // self on the reverse strand
#define SAM_FMR  32 // mate on the reverse strand
#define SAM_FR1  64 // this is read one
#define SAM_FR2 128 // this is read two
#define SAM_FSC 256 // secondary alignment

#define BWA_AVG_ERR 0.02
#define BWA_MIN_RDLEN 35 // for read trimming

#ifndef bns_pac
#define bns_pac(pac, k) ((pac)[(k)>>2] >> ((~(k)&3)<<1) & 3)
#endif

typedef int8_t	tf_flag;
typedef unsigned char ubyte_t;

typedef struct {
	char *name;
	ubyte_t *seq, *rseq;
	int len;
	int full_len;
	tf_flag rev_com; 	// Whether use reverse complement in pool
	int contig_id;
	tf_flag status;
	short pos;			// Position of some kmer on this read, used by kmer_find_reads
} bwa_seq_t;

#define BWA_MODE_GAPE       0x01
#define BWA_MODE_COMPREAD   0x02
#define BWA_MODE_LOGGAP     0x04
#define BWA_MODE_NONSTOP    0x10
#define BWA_MODE_BAM        0x20
#define BWA_MODE_BAM_SE     0x40
#define BWA_MODE_BAM_READ1  0x80
#define BWA_MODE_BAM_READ2  0x100
#define BWA_MODE_IL13       0x200

#define BWA_PET_STD   1
#define BWA_PET_SOLID 2

typedef struct {
	int max_isize, force_isize;
	int max_occ;
	int n_multi, N_multi;
	int type, is_sw, is_preload;
	double ap_prior;
} pe_opt_t;

struct __bwa_seqio_t;
typedef struct __bwa_seqio_t bwa_seqio_t;

#ifdef __cplusplus
extern "C" {
#endif

	bwa_seqio_t *bwa_seq_open(const char *fn);
	void bwa_seq_close(bwa_seqio_t *bs);
	void seq_reverse(int len, ubyte_t *seq, int is_comp);
	bwa_seq_t *bwa_read_seq(bwa_seqio_t *seq, int n_needed, int *n, int mode, int trim_qual);
	void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs);

#ifdef __cplusplus
}
#endif

#endif
