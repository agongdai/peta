/*
 * kmer.c
 *
 *  Created on: 13-Mar-2013
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <glib.h>
#include <time.h>
#include "utils.h"
#include "pehash.h"
#include "bwtaln.h"
#include "kmer.h"
#include "pehash.h"
#include "rnaseq.h"
#include "list.h"
#include "peseq.h"
#include "roadmap.h"
#include "readrm.h"
#include "mate.h"
#include "merge.h"
#include "gvdb-builder.h"

int kmer_ctg_id = 0;
int ins_size = 0;
int sd_ins_size = 0;
struct timespec kmer_start_time, kmer_finish_time;

mer *new_mer() {
	mer *m = (mer*) malloc(sizeof(mer));
	m->count = 0;
	m->s = 0;
	m->status = FRESH;
	return m;
}

mer_meta *new_mer_meta() {
	mer_meta *m = (mer_meta*) malloc(sizeof(mer_meta));
	m->n_seqs = 0;
	m->n_kmers = 0;
	m->k = 0;
	return m;
}

void free_kmer_list(GPtrArray *kmer_list) {
	int i = 0;
	mer *m = NULL;
	for (i = 0; i < kmer_list->len; i++) {
		m = g_ptr_array_index(kmer_list, i);
		free(m);
	}
	g_ptr_array_free(kmer_list, TRUE);
}

uint64_t *get_kmer_int(const ubyte_t *seq, const int start,
		const int interleaving, const int len) {
	uint64_t *key = (uint64_t*) malloc(sizeof(uint64_t));
	int i = 0;
	*key = 0;
	for (i = 0; i < len * interleaving; i += interleaving) {
		if (seq[i + start] < 0 || seq[i + start] > 4) {
			return 0;
		}
		*key <<= 2;
		// printf("%d \n", seq[i + start]);
		// Keep the last two digits and perform 'or' oepration
		*key = *key | (3 & seq[i + start]);
		// printf("%d: key = %" ID64 "\n", i, key);
	}
	return key;
}

bwa_seq_t *get_kmer_seq(mer *kmer, const int k) {
	ubyte_t *seq = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	uint64_t copy = kmer->s, index = 0;
	read = blank_seq();
	free(read->seq);
	read->seq = (ubyte_t*) calloc(k + 1, sizeof(ubyte_t));
	seq = read->seq;
	for (i = 0; i < k; i++) {
		index = copy;
		index &= 3; // Keep only the last two bits
		seq[k - 1 - i] = index;
		copy >>= 2;
	}
	read->len = k;
	read->name = (char*) malloc(32);
	sprintf(read->name, "%d", kmer->count);
	set_rev_com(read);
	return read;
}

void free_int64(gpointer data) {
	free(data);
}

void build_kmers_gvdb(const char *fa_fn, const char *out_fn, const int k) {
	bwa_seq_t *reads = NULL, *r = NULL;
	uint32_t n_reads = 0, i = 0, j = 0;
	uint64_t *mer_v = 0, n_kmers = 0;
	mer *m = NULL;
	mer_meta *meta = NULL;
	GHashTable* hash = gvdb_hash_table_new(NULL, "some");
	GPtrArray* kmer_list = g_ptr_array_sized_new(BUFSIZ);
	show_msg(__func__, "Building the kmers hashtable of %s... \n", fa_fn);
	reads = load_reads(fa_fn, &n_reads);
	for (i = 0; i < 10000; i++) {
		r = &reads[i];
		//p_query(__func__, r);
		if (i % 1000000 == 0)
			show_debug_msg(__func__,
					"%d sequences kmer counted: %d kmers...\n", i, n_kmers);
		for (j = 0; j < r->len - k; j++) {
			mer_v = get_kmer_int(r->seq, j, 1, k);
			//show_debug_msg(__func__, "Kmer id: %" ID64 " \n", *mer_v);
			m = g_hash_table_lookup(hash, mer_v);
			if (m) {
				m = (mer*) m;
				m->count++;
				free(mer_v);
				//show_debug_msg(__func__, "Kmer exists: %" ID64 ", %d \n", exist_m->s, exist_m->count);
			} else {
				m = new_mer();
				m->s = *mer_v;
				m->count++;
				n_kmers++;
				//show_debug_msg(__func__, "Kmer added: %" ID64 ", %d \n", m->s, m->count);
				g_ptr_array_add(kmer_list, m);
				g_hash_table_insert(hash, (gpointer) (mer_v), (gpointer) m);
			}
		}
	}
	meta = new_mer_meta();
	meta->k = k;
	meta->n_kmers = n_kmers;
	meta->n_seqs = n_reads;
	show_debug_msg(__func__, "%d distinct %d-mers counted \n", n_kmers, k);
	g_ptr_array_sort(kmer_list, (GCompareFunc) cmp_kmer_by_count);
	gvdb_table_write_contents(hash, out_fn, TRUE, NULL);
	free_kmer_list(kmer_list);
	g_hash_table_destroy(hash);
	bwa_free_read_seq(n_reads, reads);
}

void build_kmers(const char *fa_fn, const char *out_fn, const int k) {
	bwa_seq_t *reads = NULL, *r = NULL;
	uint32_t n_reads = 0, i = 0, j = 0;
	uint64_t *mer_v = 0, n_kmers = 0;
	mer *m = NULL;
	mer_meta *meta = NULL;
	GHashTable* hash = g_hash_table_new_full(g_int64_hash, g_int64_equal,
			(GDestroyNotify) free_int64, NULL);
	GPtrArray* kmer_list = g_ptr_array_sized_new(BUFSIZ);
	FILE *out = NULL;
	show_msg(__func__, "Building the kmers hashtable of %s... \n", fa_fn);
	reads = load_reads(fa_fn, &n_reads);
	for (i = 0; i < n_reads; i++) {
		r = &reads[i];
		//p_query(__func__, r);
		if (i % 1000000 == 0)
			show_debug_msg(__func__,
					"%d sequences kmer counted: %d kmers...\n", i, n_kmers);
		for (j = 0; j < r->len - k; j++) {
			mer_v = get_kmer_int(r->seq, j, 1, k);
			//show_debug_msg(__func__, "Kmer id: %" ID64 " \n", *mer_v);
			m = g_hash_table_lookup(hash, mer_v);
			if (m) {
				m = (mer*) m;
				m->count++;
				free(mer_v);
				//show_debug_msg(__func__, "Kmer exists: %" ID64 ", %d \n", exist_m->s, exist_m->count);
			} else {
				m = new_mer();
				m->s = *mer_v;
				m->count++;
				n_kmers++;
				//show_debug_msg(__func__, "Kmer added: %" ID64 ", %d \n", m->s, m->count);
				g_ptr_array_add(kmer_list, m);
				g_hash_table_insert(hash, (gpointer) (mer_v), (gpointer) m);
			}
		}
	}
	meta = new_mer_meta();
	meta->k = k;
	meta->n_kmers = n_kmers;
	meta->n_seqs = n_reads;
	show_debug_msg(__func__, "%d distinct %d-mers counted \n", n_kmers, k);
	g_ptr_array_sort(kmer_list, (GCompareFunc) cmp_kmer_by_count);
	out = xopen(out_fn, "w");
	fwrite(meta, sizeof(mer_meta), 1, out);
	for (i = 0; i < kmer_list->len; i++) {
		m = g_ptr_array_index(kmer_list, i);
		if (i < 10)
			show_debug_msg(__func__, "Kmer added: %" ID64 ", %d \n", m->s,
					m->count);
		fwrite(m, sizeof(mer), 1, out);
	}
	free_kmer_list(kmer_list);
	g_hash_table_destroy(hash);
	bwa_free_read_seq(n_reads, reads);
	fclose(out);
}

/**
 * Check kmer frequencies of all possible next chars, get the highest one
 * Return: ACGTN -> 01234;
 *         None -> -1;
 */
int next_char_by_kmers(GHashTable *kmers, mer_meta *meta, bwa_seq_t *query,
		const int ori) {
	bwa_seq_t *next_q = NULL;
	int counters[5], i = 0, max = 0;
	uint64_t *kmer_int = 0;
	mer *m = NULL;
	for (i = 0; i < 5; i++) {
		next_q = new_seq(query, query->len, 0);
		ext_que(next_q, i, ori);
		counters[i] = 0;
		// Check the forward kmer
		//p_query("SEQ", next_q);
		kmer_int = get_kmer_int(next_q->seq, 0, 1, meta->k);
		m = g_hash_table_lookup(kmers, kmer_int);
		if (m) {
			if (m->status != USED)
				counters[i] = m->count;
		}
		free(kmer_int);

		// Check the reverse kmer
		kmer_int = get_kmer_int(next_q->rseq, 0, 1, meta->k);
		//p_query("REV", next_q);
		if (m) {
			if (m->status != USED)
				counters[i] = m->count;
		}
		free(kmer_int);
		bwa_free_read_seq(1, next_q);
	}
	// Get max occurrence of the chars
	for (i = 0; i < 5; i++) {
		max = (counters[i] > max) ? counters[i] : max;
	}
	//show_debug_msg(__func__, "Max: %d \n", max);
	if (max == 0)
		return -1;
	for (i = 0; i < 5; i++) {
		if (max == counters[i])
			return i;
	}
	return -1;
}

/**
 * Mark a kmer as used.
 * The parameter 'kmer' is a 'copy' only, need to mark the kmer on the list
 */
void mark_kmer_used(GHashTable *kmers, mer_meta *meta, bwa_seq_t *query) {
	uint64_t *kmer_int = 0;
	mer *m = NULL;
	kmer_int = get_kmer_int(query->seq, 0, 1, meta->k);
	m = g_hash_table_lookup(kmers, kmer_int);
	if (m) {
		m->status = USED;
	}
	free(kmer_int);
	kmer_int = get_kmer_int(query->rseq, 0, 1, meta->k);
	m = g_hash_table_lookup(kmers, kmer_int);
	if (m) {
		m->status = USED;
	}
	free(kmer_int);
}

int ext_eg_query_n(GHashTable *kmers, mer_meta *meta, bwa_seq_t *query,
		edge *eg, bwa_seq_t *read, const int shift, const int ori) {
	int i = 0, ext_len = 0;
	ubyte_t *s = NULL, next_c = 0;
	if (shift > query->len)
		return 0;
	s = read->seq;
	for (i = shift; i < read->len; i++) {
		next_c = s[i];
		ext_que(query, next_c, ori);
		ext_con(eg->contig, next_c, 0);
		//p_ctg_seq(__func__, eg->contig);
		mark_kmer_used(kmers, meta, query);
	}
	ext_len = eg->contig->len - eg->len;
	eg->len = eg->contig->len;
	return ext_len;
}

int ext_by_mates(GHashTable *kmers, mer_meta *meta, edgearray *mates,
		bwa_seq_t *query, edge *eg, const int ori) {
	int ol = 0, i = 0, ext_len = 0;
	bwa_seq_t *m = NULL, *tmp = NULL, *template = NULL;
	template = new_seq(eg->contig, query->len, eg->len - query->len);
	if (ori) {
		seq_reverse(template->len, template->seq, 0);
	}
	for (i = 0; i < mates->len; i++) {
		m = g_ptr_array_index(mates, i);
		tmp = m;
		if (m->rev_com)
			tmp = new_mem_rev_seq(m, m->len, 0);
		ol = find_ol_within_k(tmp, template, 1, 0, template->len, ori);
		//p_ctg_seq(__func__, tmp);
		//p_ctg_seq(__func__, template);
		//show_debug_msg(__func__, "OL: %d \n", ol);
		if (ol >= 8) {
			//p_query("Mate", m);
			ext_len = ext_eg_query_n(kmers, meta, query, eg, tmp, ol, ori);
		}

		if (m->rev_com)
			bwa_free_read_seq(1, tmp);
		if (ext_len > 0)
			break;
	}
	bwa_free_read_seq(1, template);
	return ext_len;
}

void kmer_ext_edge(edge *eg, bwa_seq_t *query, GHashTable *kmers,
		mer_meta *meta, hash_table *ht, const int ori) {
	int next_c = 0, pre_round_len = eg->len;
	pool *mate_pool = NULL;
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	p_query(__func__, query);
	while (1) {
		//p_query(__func__, query);
		//p_ctg_seq("Contig", eg->contig);
		next_c = next_char_by_kmers(kmers, meta, query, ori);
		//show_debug_msg(__func__, "Next char: %d \n", next_c);
		if (next_c == -1) {
			//show_debug_msg(__func__, "Realining edge [%d, %d] \n", eg->id,
			//		eg->len);
			//p_ctg_seq("Contig", eg->contig);
			realign_extended(ht, eg, pre_round_len, 0, ori);
			//p_readarray(eg->reads, 1);
			//show_debug_msg(__func__, "Getting mate pool \n");
			mate_pool = get_mate_pool_from_edge(eg, ht, ori, ins_size,
					sd_ins_size);
			//p_readarray(mate_pool->reads, 1);
			show_debug_msg(__func__, "Extending by mates \n");
			if (ext_by_mates(kmers, meta, mate_pool->reads, query, eg, ori)) {
				pre_round_len = eg->len;
				continue;
			} else
				break;
		}
		ext_con(eg->contig, next_c, 0);
		eg->len = eg->contig->len;
		ext_que(query, next_c, ori);
		mark_kmer_used(kmers, meta, query);
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Edge [%d, %d] \n", eg->id, eg->len);
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	bwa_free_read_seq(1, query);
}

/**
 * Extend a kmer
 */
edge *kmer_ext(GHashTable *kmers, mer_meta *meta, bwa_seq_t *kmer,
		hash_table *ht) {
	edge *eg = NULL;
	bwa_seq_t *query = NULL;
	int round_1_len = 0, round_2_len = 0;

	eg = new_eg();
	eg->id = kmer_ctg_id++;
	// Get a copy of the kmer
	query = new_seq(kmer, kmer->len, 0);
	eg->contig = new_seq(query, query->len, 0);
	kmer_ext_edge(eg, query, kmers, meta, ht, 0);
	round_1_len = eg->len;
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	query = new_seq(eg->contig, kmer->len, 0);
	kmer_ext_edge(eg, query, kmers, meta, ht, 1);
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		query = new_seq(eg->contig, kmer->len, 0);
		round_2_len = eg->len;
		kmer_ext_edge(eg, query, kmers, meta, ht, 0);
		show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			query = new_seq(eg->contig, kmer->len, 0);
			kmer_ext_edge(eg, query, kmers, meta, ht, 1);
			show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		}
	}

	return eg;
}

GHashTable *load_kmers(const char *kmer_file, GPtrArray *kmer_list,
		mer_meta *meta) {
	FILE *kmers_f = xopen(kmer_file, "r");
	uint64_t i = 0, reading = 0;
	mer *m = NULL;
	GHashTable* hash = g_hash_table_new_full(g_int64_hash, g_int64_equal,
			(GDestroyNotify) free_int64, NULL);
	show_msg(__func__, "Loading kmers...\n");
	reading = fread(meta, sizeof(mer_meta), 1, kmers_f);
	show_msg(__func__, "%" ID64 " kmers to load... \n", meta->n_kmers);
	for (i = 0; i < meta->n_kmers; i++) {
		m = new_mer();
		reading = fread(m, sizeof(mer), 1, kmers_f);
		g_hash_table_insert(hash, (gpointer) (&(m->s)), (gpointer) m);
		if (m->count > 2)
			g_ptr_array_add(kmer_list, m);
	}
	fclose(kmers_f);
	return hash;
}

void ext_by_kmers(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads) {
	hash_table *ht = NULL;
	uint64_t i = 0;
	GPtrArray *all_edges = NULL, *kmer_list = NULL;
	mer *m = NULL;
	mer_meta *meta = NULL;
	bwa_seq_t *kmer_seq = NULL;
	edge *eg = NULL;
	FILE *kmer_contigs = NULL;
	GHashTable *kmers = NULL;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	ins_size = insert_size;
	sd_ins_size = sd_insert_size;

	kmer_list = g_ptr_array_sized_new(BUFSIZ);
	meta = (mer_meta*) malloc(sizeof(mer_meta));
	kmers = load_kmers(kmer_file, kmer_list, meta);
	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	all_edges = g_ptr_array_sized_new(BUFSIZ); // kmer_list->kmers->len
	ht = pe_load_hash(lib_file);
	show_msg(__func__, "Extending by kmers...\n");
	for (i = 0; i < kmer_list->len; i++) {
		m = g_ptr_array_index(kmer_list, i);
		if (m->status != USED && m->count > 2) {
			kmer_seq = get_kmer_seq(m, meta->k);
			p_query(__func__, kmer_seq);
			if (has_n(kmer_seq, 4) || is_biased_q(kmer_seq)) {
				bwa_free_read_seq(1, kmer_seq);
				continue;
			}
			//if (all_edges->len > 100)
			//	break;
			show_debug_msg(__func__,
					"========================== %d ===================== \n", i);
			eg = kmer_ext(kmers, meta, kmer_seq, ht);
			if (eg->len > 100)
				g_ptr_array_add(all_edges, eg);
			bwa_free_read_seq(1, kmer_seq);
			//break;
		}
	}
	kmer_contigs = xopen("../SRR097897_out/kmer_contigs.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
	g_ptr_array_free(kmer_list, TRUE);

	reset_edge_ids(all_edges);
	merge_ol_edges(all_edges, insert_size, sd_insert_size, ht, n_threads);
	destroy_ht(ht);
	kmer_contigs = xopen("../SRR097897_out/merged.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
}
