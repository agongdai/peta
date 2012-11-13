/*
 * pealn.c
 *
 *  Created on: 26-Jun-2011
 *      Author: carl
 */

#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "pealn.h"
#include "pehash.h"
#include "bwase.h"
#include "pechar.h"

void free_hits(GPtrArray *hits) {
	pos_tuple *h;
	int i = 0;
	if (hits && hits->len > 0) {
		for (i = 0; i < hits->len; i++) {
			h = g_ptr_array_index(hits, i);
			free(h);
		}
	}
	g_ptr_array_free(hits, TRUE);
}

int compare_hits(gpointer a, gpointer b) {
	pos_tuple *h1 = *((pos_tuple**) a);
	pos_tuple *h2 = *((pos_tuple**) b);
	return h1->seq_id - h2->seq_id;
}

void quicksort(int *arr, int elements) {
#define  MAX_LEVELS  300
	int piv, beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, left, right, swap;
	beg[0] = 0;
	end[0] = elements;
	while (i >= 0) {
		left = beg[i];
		right = end[i] - 1;
		if (left < right) {
			piv = arr[left];
			while (left < right) {
				while (arr[right] >= piv && left < right)
					right--;
				if (left < right)
					arr[left++] = arr[right];
				while (arr[left] <= piv && left < right)
					left++;
				if (left < right)
					arr[right--] = arr[left];
			}
			arr[left] = piv;
			beg[i + 1] = left + 1;
			end[i + 1] = end[i];
			end[i++] = left;
			if (end[i] - beg[i] > end[i - 1] - beg[i - 1]) {
				swap = beg[i];
				beg[i] = beg[i - 1];
				beg[i - 1] = swap;
				swap = end[i];
				end[i] = end[i - 1];
				end[i - 1] = swap;
			}
		} else {
			i--;
		}
	}
}

void free_alg(alignarray *alns) {
	int i = 0;
	alg *a;
	if (alns) {
		for (i = 0; i < alns->len; i++) {
			a = g_ptr_array_index(alns, i);
			free(a);
		}
		g_ptr_array_free(alns, TRUE);
	}
}

void reset_alg(alignarray *alns) {
	int i = 0;
	alg *a;
	if (alns) {
		for (i = 0; i < alns->len; i++) {
			a = g_ptr_array_index(alns, i);
			// printf("[reset_alg] %p: %" ID64 ": %d, %" ID64 "\n", a, a->r_id, a->pos, a->q_id);
			free(a);
		}
		while (alns->len > 0)
			g_ptr_array_remove_index(alns, 0);
	}
}

alg *aligned(const alignarray *alns, const index64 id) {
	unsigned int i = 0;
	alg *a;
	for (i = 0; i < alns->len; i++) {
		a = g_ptr_array_index(alns, i);
		if (a->r_id == id)
			return a;
	}
	return 0;
}

void p_align(const alignarray *alns) {
	unsigned int i = 0;
	alg *a;
	printf("[p_align]****************************** # of alignments: %d\n",
			alns->len);
	for (i = 0; i < alns->len; i++) {
		a = g_ptr_array_index(alns, i);
		printf("[p_align] %d: %" ID64 ", %d \n", i, a->r_id, a->pos);
	}
	printf("[p_align]****************************** \n");
}

pos_tuple *set_hit(const index64 seq_id, const int shift, const int offset) {
	pos_tuple *h = (pos_tuple*) malloc(sizeof(pos_tuple));
	h->seq_id = seq_id;
	h->offset = offset;
	h->shift = shift;
	return h;
}

alg *new_alg() {
	alg *a = (alg*) malloc(sizeof(alg));
	a->diff = 0;
	a->flag = 0;
	a->len = 0;
	a->pos = 0;
	a->q_id = 0;
	a->r_id = 0;
	a->strand = 0;
	a->rev_comp = 0;
	return a;
}

void extract_alns(const bwa_seq_t *seqs, const bwa_seq_t *query,
		const ubyte_t *q_seq, GPtrArray *hits, const int mismatches,
		const int ol, const int is_rc, alignarray *aligns) {
	pos_tuple *h, *h2;
	alg *a;
	int i = 0, diff = 0;
	bwa_seq_t *aln_read;
	// p_query(query);
	for (i = 0; i < hits->len; i++) {
		h = g_ptr_array_index(hits, i);
		if (i > 0)
			h2 = g_ptr_array_index(hits, i - 1);
		//		printf("Alignment %d: (%" ID64 ", %d, %d) \n", i, h->seq_id, h->shift,
		//				h->offset);
		if ((i > 0 && h->seq_id == h2->seq_id)) {
			continue;
		} else {
			aln_read = &seqs[h->seq_id];
			diff = is_sub_seq_aln(q_seq, query->len, h->shift, h->offset,
					aln_read, mismatches, ol);
			if (diff != NOT_FOUND) {
				a = new_alg();
				a->diff = diff;
				a->q_id = get_index(query->name);
				a->r_id = get_index(aln_read->name);
				a->len = query->len;
				a->pos = h->shift;
				a->rev_comp = is_rc;
				g_ptr_array_add(aligns, a);
			}
		}
	}
}

void pe_aln_query(const bwa_seq_t *query, const ubyte_t *q_seq,
		const hash_table *ht, const int mismatches, const int ol,
		const int is_rc, alignarray *aligns) {
	index64 seq_id = 0, pos_index = 0;
	int i = 0, n_occ = 0, j = 0, locus = 0, n_hits = 0;
	hash_opt *hash_o;
	hash_key *k_mers_occ_acc, key;
	hash_value *pos, value;
	GPtrArray *hits;
	pos_tuple *h;
	bwa_seq_t *seqs;

	seqs = ht->seqs;
	hash_o = ht->o;
	k_mers_occ_acc = ht->k_mers_occ_acc;
	pos = ht->pos;
	hits = g_ptr_array_sized_new(N_DEFAULT_HITS);
	for (i = 0; i <= (query->len - 2 * hash_o->k + 1); i++) {
		key = get_hash_key(q_seq, i, hash_o->interleaving, hash_o->k);
		// If contains 'n', just ignore it.
		if (!key)
			continue;
		n_occ = k_mers_occ_acc[key + 1] - k_mers_occ_acc[key];
		pos_index = k_mers_occ_acc[key];
		for (j = 0; j < n_occ; j++) {
			value = pos[pos_index + j];
			//			printf("i: %d; pos: %" ID64 "\n", i, value);
			read_hash_value(&seq_id, &locus, value);
			//			printf("id: %" ID64 ", locus: %d \n", seq_id, locus);
			h = set_hit(seq_id, locus - i, locus);
//			if (query->len == ol) { // If it is a full-length query, only return hits at pos 0
//				if (h->shift == 0) {
//					g_ptr_array_add(hits, h);
//					n_hits++;
//				} else {
//					free(h);
//				}
//			} else {
				g_ptr_array_add(hits, h);
				n_hits++;
//			}
			// printf("[pe_aln_query] K-mer starts at %d: %" ID64 ", %d \n", i,
			//		seq_id, locus);
		}
	}

	// printf("After sorting: ---------------------------------- \n");
	g_ptr_array_sort(hits, (GCompareFunc) compare_hits);
	extract_alns(seqs, query, q_seq, hits, mismatches, ol, is_rc, aligns);
	free_hits(hits);
}

void pe_aln_dummy(const char *fa_fn, const char *query_fn, const int mismatches) {
	bwa_seq_t *queries = 0, *q;
	bwa_seqio_t *ks;
	index64 n_qs = 0;
	hash_table *ht;
	clock_t t = clock();
	alignarray *alns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	int i = 0;

	ht = pe_load_hash(fa_fn);
	ks = bwa_open_reads(BWA_MODE, query_fn);
	while ((queries = bwa_read_seq(ks, 0x400000, &n_qs, BWA_MODE, 0)) != 0) {
		pe_reverse_seqs(queries, n_qs);
		break;
	}

	for (i = 0; i < n_qs; i++) {
		printf(
				"----------------------------------------------------------------\n");
		q = &queries[i];
		p_query(__func__, q);
		pe_aln_query(q, q->seq, ht, mismatches, 21, 0, alns);
		p_align(alns);
		p_query(__func__, q);
		pe_aln_query(q, q->seq, ht, mismatches, 25, 0, alns);
		p_align(alns);
	}
	bwa_seq_close(ks);
	fprintf(stderr, "[pe_aln_dummy] Done. %" ID64 " query. %.2f sec. \n", n_qs,
			(float) (clock() - t) / CLOCKS_PER_SEC);
}

int pe_aln_test(int argc, char *argv[]) {
	int c;
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 2 > argc) {
		return 1;
	}
	pe_aln_dummy(argv[optind], argv[optind + 1], 2);
	return 0;
}
