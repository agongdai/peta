#ifndef RNASEQ_H_
#define RNASEQ_H_
#include "bwase.h"
#include <inttypes.h>

#define N_CHUNK_SEQS	4194304 // # of reads read in every time
#define N_DF_MAX_SEQS	8388608
#define BWA_MODE		3

bwa_seq_t *load_reads(const char *fa_fn, uint32_t *n_seqs);

#endif /* RNASEQ_H_ */
