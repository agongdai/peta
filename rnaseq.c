#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include "bwase.h"
#include "peseq.h"
#include "utils.h"
#include "clean.h"
#include "pehash.h"

bwa_seq_t *load_reads(char *fa_fn, uint32_t *n_seqs) {
	bwa_seq_t *seqs, *part_seqs;
	bwa_seqio_t *ks;
	uint32_t n_part_seqs = 0, n_seqs_full = 0, n_seqs_loaded = 0;
	clock_t t = clock();

	ks = bwa_open_reads(BWA_MODE, fa_fn);
	n_seqs_full = N_CHUNK_SEQS;
	show_msg(__func__, "Loading reads from library %s...\n", fa_fn);
	seqs = (bwa_seq_t*) calloc (N_DF_MAX_SEQS, sizeof(bwa_seq_t));
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, BWA_MODE, 0))
			!= 0) {
		show_msg(__func__, "%d sequences loaded: %.2f sec... \n",
				n_seqs_loaded + n_part_seqs, fa_fn, (float) (clock() - t) / CLOCKS_PER_SEC);
		pe_reverse_seqs(part_seqs, n_part_seqs);

		if ((n_seqs_loaded + n_part_seqs) > n_seqs_full) {
			n_seqs_full += n_part_seqs + 2;
			kroundup32(n_seqs_full);
			seqs = (bwa_seq_t*) realloc(seqs, sizeof(bwa_seq_t) * n_seqs_full);
		}
		memmove(&seqs[n_seqs_loaded], part_seqs, sizeof(bwa_seq_t) * n_part_seqs);
		free(part_seqs);
		n_seqs_loaded += n_part_seqs;
	}
	bwa_seq_close(ks);
	if (n_seqs_loaded < 1) {
		err_fatal(__func__,
				"No sequence in file %s, make sure the format is correct! \n",
				fa_fn);
	}
	*n_seqs = n_seqs_loaded;
	return seqs;
}
