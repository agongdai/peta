/*
 * hash.c
 *
 *  Created on: 09-May-2013
 *      Author: carl
 */
#include <vector>
#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <ctime>
#include <stdlib.h>
#include <stdio.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "hash.hpp"
#include "utils.h"
#include "rnaseq.h"
#include "tpl.hpp"
#include "peseq.h"
#include "k_hash.h"

/**
 * Return template id
 */
int read_tpl_occ(uint64_t value, int *locus, int *rev_com, int *head_or_tail) {
	int t_id = value >> 32;
	*locus = value & TPL_LOCUS_LOWER;
	*rev_com = value >> 30 & 3;
	*head_or_tail = value >> 28 & 3;
	return t_id;
}

gint cmp_kmers_by_count(gpointer a, gpointer b) {
	kmer_counter *c_a = *((kmer_counter**) a);
	kmer_counter *c_b = *((kmer_counter**) b);
	return ((c_b->count) - c_a->count);
}

/**
 * 64 bits:
 * First 32: 	template id
 * Next 2:		rev_com
 * Next 2:		head_or_tail
 * Last 28:		locus
 */
uint64_t hash_tpl_occ(tpl *t, int locus, int rev_com, int head_or_tail) {
	uint64_t occ = 0;
	occ += t->id;
	occ <<= 32;
	occ += rev_com << 30;
	occ += head_or_tail << 28;
	occ += locus;
	//show_debug_msg(__func__,
	//		"id: %d; locus: %d; rev_com: %d; head_or_tail: %d => %" ID64 "\n",
	//		t->id, locus, rev_com, head_or_tail, occ);
	return occ;
}

/**
 * Build the hash table for head/tail of templates
 */
void build_tpl_hash(kmer_hash &hash, tpl_hash *tpls, const int k,
		const int read_len) {
	uint32_t i = 0, j = 0, start = 0, end = 0;
	uint64_t kmer = 0, *kmer_freq = NULL, *occ = NULL, value = 0;
	int id = 0, locus = 0, rev_com = 0, head_or_tail = 0;
	tpl *t = NULL;
	mer_counter::iterator it;
	tpl_hash::iterator im;
	mer_counter counter;
	// Count the k-mer frequencies at the head/tail of all templates
	for (im = tpls->begin(); im != tpls->end(); ++im) {
		t = (tpl*) im->second;
		p_tpl(t);
		if (t->len < k)
			continue;
		// The head
		end = t->len > read_len ? (read_len - k) : (t->len - k);
		for (j = 0; j <= end; j++) {
			kmer = get_kmer_int(t->ctg->seq, j, 1, k);
			counter[kmer]++;
			kmer = get_kmer_int(t->ctg->rseq, j, 1, k);
			counter[kmer]++;
		}
		// The tail
		start = t->len > read_len ? (t->len - read_len) : 0;
		for (j = start; j <= t->len - k; j++) {
			kmer = get_kmer_int(t->ctg->seq, j, 1, k);
			counter[kmer]++;
			kmer = get_kmer_int(t->ctg->rseq, j, 1, k);
			counter[kmer]++;
		}
	}
	for (it = counter.begin(); it != counter.end(); ++it) {
		kmer_freq = (uint64_t*) calloc(it->second + 1, sizeof(uint64_t));
		hash[it->first] = kmer_freq;
		//		show_debug_msg(__func__, "kmer: %" ID64 " -> %" ID64 ". \n", it->first,
		//				it->second);
	}
	counter.clear();
	for (im = tpls->begin(); im != tpls->end(); ++im) {
		t = (tpl*) im->second;
		if (t->len < k)
			continue;
		// The head
		end = t->len > read_len ? (read_len - k) : (t->len - k);
		for (j = 0; j <= end; j++) {
			kmer = get_kmer_int(t->ctg->seq, j, 1, k);
			kmer_freq = hash[kmer];
			kmer_freq[++kmer_freq[0]] = hash_tpl_occ(t, j, 0, 0);
			kmer = get_kmer_int(t->ctg->rseq, j, 1, k);
			kmer_freq = hash[kmer];
			kmer_freq[++kmer_freq[0]] = hash_tpl_occ(t, j, 1, 0);
		}
		// The tail
		start = t->len > read_len ? (t->len - read_len) : 0;
		for (j = start; j <= t->len - k; j++) {
			kmer = get_kmer_int(t->ctg->seq, j, 1, k);
			kmer_freq = hash[kmer];
			kmer_freq[++kmer_freq[0]] = hash_tpl_occ(t, j, 0, 1);
			kmer = get_kmer_int(t->ctg->rseq, j, 1, k);
			kmer_freq = hash[kmer];
			kmer_freq[++kmer_freq[0]] = hash_tpl_occ(t, j, 1, 1);
		}
	}
//	for (kmer_hash::iterator m = hash.begin(); m != hash.end(); ++m) {
//		kmer = m->first;
//		occ = (uint64_t*) m->second;
//		if (occ[0] > 0) {
//			bwa_seq_t *key_seq = get_key_seq(kmer, 11);
//			p_query(__func__, key_seq);
//			bwa_free_read_seq(1, key_seq);
//			for (i = 0; i < occ[0]; i++) {
//				value = occ[i + 1];
//				show_debug_msg(__func__, "VALUE: %" ID64 "\n", value);
//				id = read_tpl_occ(value, &locus, &rev_com, &head_or_tail);
//				show_debug_msg(__func__,
//						"Template %d, locus %d, rev_com %d, head_or_tail %d \n --- \n",
//						id, locus, rev_com, head_or_tail);
//			}
//		}
//	}
}
