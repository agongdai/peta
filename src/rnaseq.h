#ifndef RNASEQ_H_
#define RNASEQ_H_
#include "bwtaln.h"
#include <inttypes.h>

#define N_CHUNK_SEQS	4194304 // # of reads read in every time
#define N_DF_MAX_SEQS	8388608
#define BWA_MODE		3

extern unsigned char nst_nt4_table[256];

bwa_seq_t *load_reads(const char *fa_fn, uint32_t *n_seqs);
bwa_seq_t *load_arr_reads(const char *fa_fn, uint32_t *n_reads);

#endif /* RNASEQ_H_ */
