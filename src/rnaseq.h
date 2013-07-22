#ifndef RNASEQ_H_
#define RNASEQ_H_
#include "bwtaln.h"
#include "k_hash.h"
#include <glib.h>
#include <inttypes.h>


#define N_DF_MAX_SEQS	8388608

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char nst_nt4_table[256];

bwa_seq_t *load_reads(const char *fa_fn, uint64_t *n_seqs);
bwa_seq_t *load_arr_reads(const char *fa_fn, uint32_t *n_reads);
GPtrArray *load_solid_reads(const char *solid_fn, bwa_seq_t *seqs,
		const int n_seqs);

#ifdef __cplusplus
}
#endif

#endif /* RNASEQ_H_ */
