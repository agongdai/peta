#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include "peseq.h"
#include "rnaseq.h"

// Encode the sequence into 2-bit representations:
// 		a: 00; c: 01; g: 10; t: 11
// In the ascii table, a/A is 97/65; c/C is 99/67; g/G is 103/71; t/T is 116/84
unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

bwa_seq_t *load_reads(const char *fa_fn, uint32_t *n_reads) {
	bwa_seq_t *seqs, *part_seqs;
	bwa_seqio_t *ks;
	int n_part_seqs = 0, n_seqs_full = 0, n_seqs_loaded = 0;
	clock_t t = clock();

	ks = bwa_open_reads(BWA_MODE, fa_fn);
	n_seqs_full = N_CHUNK_SEQS;
	show_msg(__func__, "Loading reads from library %s...\n", fa_fn);
	seqs = (bwa_seq_t*) calloc(N_DF_MAX_SEQS, sizeof(bwa_seq_t));
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, BWA_MODE,
			0)) != 0) {
		show_msg(__func__, "%d sequences loaded: %.2f sec... \n", n_seqs_loaded
				+ n_part_seqs, fa_fn, (float) (clock() - t) / CLOCKS_PER_SEC);
		pe_reverse_seqs(part_seqs, n_part_seqs);
		if ((n_seqs_loaded + n_part_seqs) > n_seqs_full) {
			n_seqs_full += n_part_seqs + 2;
			kroundup32(n_seqs_full);
			seqs = (bwa_seq_t*) realloc(seqs, sizeof(bwa_seq_t) * n_seqs_full);
		}
		memmove(&seqs[n_seqs_loaded], part_seqs, sizeof(bwa_seq_t)
				* n_part_seqs);
		free(part_seqs);
		n_seqs_loaded += n_part_seqs;
	}
	bwa_seq_close(ks);
	if (n_seqs_loaded < 1) {
		err_fatal(__func__,
				"No sequence in file %s, make sure the format is correct! \n",
				fa_fn);
	}
	*n_reads = n_seqs_loaded;
	return seqs;
}

/**
 * Load reads in a raw way
 */
bwa_seq_t *load_arr_reads(const char *fa_fn, uint32_t *n_reads) {
	bwa_seq_t *reads = NULL, *r = NULL;
	int h_len = 0, len = 0, len_full = 0;
	uint32_t i = 0, n_full = 0, n_loaded = 0;
	clock_t t = clock();
	FILE *fa = NULL;
	fa = xopen(fa_fn, "r");
	char c = 0;
	short in_header = 0, reading_h = 0;
	char *header = (char*) calloc(BUFSIZ, sizeof(char));
	ubyte_t *seq = NULL;

	reads = (bwa_seq_t*) calloc(N_CHUNK_SEQS, sizeof(bwa_seq_t));
	if (!reads)
		err_fatal(__func__, "Failed to allocate memory! Exit. \n");
	while ((c = fgetc(fa)) != EOF) {
		if (c == '>') {
			if (n_loaded % N_CHUNK_SEQS == 0) {
				show_msg(__func__, "%d reads loaded, %d full...\n", n_loaded,
						n_full);
			}
			in_header = 1;
			reading_h = 1;
			if ((n_loaded + 2) > n_full) {
				n_full += 2;
				kroundup32(n_full);
				reads = (bwa_seq_t*) realloc(reads, sizeof(bwa_seq_t) * n_full);
				if (!reads) {
					err_fatal(__func__, "Failed to allocate memory! Exit. \n");
				}
			}
			// The value of 'reads' may be changed to another pointer after realloc!
			r = &reads[n_loaded++];
			h_len = 0;
			len_full = 60;
			len = 0;
			seq = (ubyte_t*) calloc(len_full, sizeof(ubyte_t));
			r->seq = seq;
			r->rseq = NULL;
			continue;
		}
		if (c == '\n') {
			if (in_header) {
				header[h_len] = '\0';
				r->name = strdup(header);
				in_header = 0;
			}
			continue;
		}
		if (in_header) {
			if (c == ' ')
				reading_h = 0;
			if (reading_h)
				header[h_len++] = c;
		} else {
			seq[len++] = (ubyte_t) nst_nt4_table[(c)];
			r->len = len;
			if ((len + 2) > len_full) {
				len_full += 2;
				kroundup32(len_full);
				seq = (ubyte_t*) realloc(seq, sizeof(ubyte_t) * len_full);
				if (!seq) {
					err_fatal(__func__, "Failed to allocate memory! Exit. \n");
				}
				r->seq = seq;
			}
		}
	}
	reads = (bwa_seq_t*) realloc(reads, n_loaded * sizeof(bwa_seq_t));
	*n_reads = n_loaded;
	for (i = 0; i < n_loaded; i++) {
		r = &reads[i];
		set_rev_com(r);
	}

	show_msg(__func__, "%d sequences loaded in %.2f sec...", n_loaded,
			(float) (clock() - t) / CLOCKS_PER_SEC);

	free(header);
	fclose(fa);
	return reads;
}

GPtrArray *load_solid_reads(const char *solid_fn, bwa_seq_t *seqs,
		const int n_seqs) {
	int i = 0;
	char line[80];
	GPtrArray *solid_reads = NULL;
	bwa_seq_t *query = NULL;
	FILE *solid = xopen(solid_fn, "r");

	solid_reads = g_ptr_array_sized_new(n_seqs / 10);
	while (fgets(line, 80, solid) != NULL) {
		i = atoi(line);
		query = &seqs[i];
		g_ptr_array_add(solid_reads, query);
	}
	fclose(solid);
	return solid_reads;
}
