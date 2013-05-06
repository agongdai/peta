/*
 * ass.cpp
 *
 *  Created on: 09-Apr-2013
 *      Author: carl
 */

#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "edge.h"
#include "edgelist.h"
#include "mate.h"
#include "pool.h"
#include "merge.h"
#include "readrm.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "pelib.h"
#include "utils.h"
#include "pehash.h"
#include "peseq.h"
#include "clean.h"

using namespace std;

int kmer_ctg_id = 1;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;
GPtrArray *branching_events = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

edge *blank_edge(uint64_t query_int, int kmer_len, int init_len, int ori) {
	bwa_seq_t *kmer = NULL;
	edge *eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	kmer = get_kmer_seq(query_int, kmer_len);
	p_query(__func__, kmer);
	if (ori)
		eg->contig = new_seq(kmer, init_len, 0);
	else
		eg->contig = new_seq(kmer, init_len, kmer->len - init_len);
	eg->len = eg->contig->len;
	eg->tid = atoi(kmer->name);
	eg->name = strdup(kmer->name);
	bwa_free_read_seq(1, kmer);
	return eg;
}

void add_a_junction(edge *main_tpl, edge *branch_tpl, uint64_t kmer, int locus,
		int ori, int weight) {
	junction *j = (junction*) malloc(sizeof(junction));
	j->main_tpl = main_tpl;
	j->branch_tpl = branch_tpl;
	j->locus = main_tpl->len;
	j->ori = ori;
	j->locus = locus;
	j->weight = weight;
	j->kmer = kmer;
	g_mutex_lock(kmer_id_mutex);
	g_ptr_array_add(branching_events, j);
	g_mutex_unlock(kmer_id_mutex);
}

GPtrArray *hash_to_array(tpl_hash *all_tpls) {
	uint64_t id = 0;
	edge *eg = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		id = m->first;
		eg = (edge*) m->second;
		g_ptr_array_add(tpls, eg);
	}
	return tpls;
}

gint cmp_kmers_by_count(gpointer a, gpointer b) {
	kmer_counter *c_a = *((kmer_counter**) a);
	kmer_counter *c_b = *((kmer_counter**) b);
	return ((c_b->count) - c_a->count);
}

void mark_tpl_kmers_used(edge *eg, hash_map *hm, const int kmer_len) {
	int i = 0;
	uint64_t query_int = 0;
	for (i = 0; i < eg->len - kmer_len; i++) {
		query_int = get_kmer_int(eg->contig->seq, i, 1, kmer_len);
		mark_kmer_used(query_int, hm, eg->id, i);
	}
}

void mark_tpl_kmers_fresh(edge *eg, hash_map *hm, const int kmer_len) {
	int i = 0;
	uint64_t query_int = 0;
	for (i = 0; i < eg->len - kmer_len; i++) {
		query_int = get_kmer_int(eg->contig->seq, i, 1, kmer_len);
		mark_kmer_not_used(query_int, hm);
	}
}

/**
 * Update the junction locus for those edges right connected to itself.
 * Because the locus is not correct when the junction is recorded.
 */
void upd_tpl_jun_locus(edge *eg, GPtrArray *branching_events,
		const uint32_t kmer_len) {
	uint32_t i = 0, k = 0;
	uint64_t query_int = 0;
	junction *jun = NULL;
	for (k = 0; k < branching_events->len; k++) {
		jun = (junction*) g_ptr_array_index(branching_events, k);
		if (jun->main_tpl == eg && jun->branch_tpl == eg) {
			for (i = 0; i < eg->len - kmer_len; i++) {
				query_int = get_kmer_int(eg->contig->seq, i, 1, kmer_len);
				if (query_int == jun->kmer) {
					jun->locus = i;
					break;
				}
			}
		}
	}
}

/**
 * Extend the short template as int64
 * max_len: should be smaller than 64
 */
uint64_t try_short_tpl_ext(hash_map *hm, uint64_t query, const int ori,
		const int max_len) {
	uint64_t branch_template = 0;
	int i = 0, c = 0;
	for (i = 0; i < max_len; i++) {
		c = next_char_by_kmers(hm, query, 1, ori);
		if (c == -1)
			return LOWER_ONES_32;
		else {
			branch_template = shift_bit(branch_template, c, hm->o->k, ori);
			query = shift_bit(query, c, hm->o->k, ori);
			//show_debug_msg(__func__, "Branch int: %"ID64"\n", branch_template);
		}
	}
	return branch_template;
}

/**
 * Extend the short template as bwa_seq_t, max_len may be larger than 64
 */
bwa_seq_t *try_short_tpl_byte_ext(hash_map *hm, uint64_t query, int first_c,
		const int ori, const int max_len) {
	bwa_seq_t *branch_seq = NULL;
	int i = 0, c = 0;
	branch_seq = blank_seq(max_len);
	branch_seq->seq[0] = first_c;
	branch_seq->len = 1;
	for (i = 1; i < max_len; i++) {
		c = next_char_by_kmers(hm, query, 1, ori);
		if (c == -1)
			break;
		branch_seq->seq[i] = c;
		branch_seq->len = i + 1;
		query = shift_bit(query, c, hm->o->k, ori);
		//p_ctg_seq(__func__, branch_seq);
	}
	if (ori)
		seq_reverse(branch_seq->len, branch_seq->seq, 0);
	set_rev_com(branch_seq);
	return branch_seq;
}

/**
 * Check whether there are reads in the junction area
 */
int find_junction_reads(hash_map *hm, bwa_seq_t *left, bwa_seq_t *right,
		const int max_len, int *weight) {
	int left_len = 0, right_len = 0, n_reads = 0;
	GPtrArray *reads = NULL;
	bwa_seq_t *junction_seq = blank_seq(max_len);

	left_len = (left->len > max_len / 2) ? (max_len / 2) : left->len;
	memcpy(junction_seq->seq, left->seq + (left->len - left_len),
			sizeof(ubyte_t) * left_len);
	right_len = (right->len) > (max_len / 2) ? (max_len / 2) : (right->len);
	memcpy(junction_seq->seq + left_len, right->seq, sizeof(ubyte_t)
			* right_len);
	junction_seq->len = left_len + right_len;
	set_rev_com(junction_seq);
	p_query("Junction seq", junction_seq);
	reads = kmer_find_reads(junction_seq, hm, 1, 0);
	n_reads = reads->len;
	show_debug_msg(__func__, "# of junction reads: %d \n", n_reads);
	//p_readarray(reads, 1);
	*weight = n_reads;
	bwa_free_read_seq(1, junction_seq);
	g_ptr_array_free(reads, TRUE);
	if (n_reads > 0)
		return 1;
	return 0;
}

/**
 * If there are reads spanning the junction, the junction is valid
 * Parameter is_rev: the input edge is reversed
 */
int val_junction_reads(hash_map *hm, edge *main_tpl, uint64_t query_int, int c,
		const int ori, int *weight, const int max_len) {
	uint64_t query = 0, is_valid = 1;
	bwa_seq_t *branch_seq = NULL, *main_tail = NULL;
	query = shift_bit(query_int, c, hm->o->k, ori);
	branch_seq = try_short_tpl_byte_ext(hm, query, c, ori, max_len / 2);
	if (ori) {
		seq_reverse(main_tpl->len, main_tpl->contig->seq, 0);
		main_tail = cut_edge_tail(main_tpl, 0, max_len / 2, ori);
	} else
		main_tail = cut_edge_tail(main_tpl, main_tpl->len, max_len / 2, ori);
	p_ctg_seq("Main", main_tail);
	p_ctg_seq("Branch", branch_seq);
	if (ori)
		is_valid = find_junction_reads(hm, branch_seq, main_tail, max_len,
				weight);
	else
		is_valid = find_junction_reads(hm, main_tail, branch_seq, max_len,
				weight);
	bwa_free_read_seq(1, branch_seq);
	bwa_free_read_seq(1, main_tail);
	if (ori)
		seq_reverse(main_tpl->len, main_tpl->contig->seq, 0);
	return is_valid;
}

/**
 Validate the branching event on 'main_tpl' at locus 'shift' with orientation 'ori'
 The parameter c is the first char of the candidate branch
 **/
int val_branching(hash_map *hm, edge *main_tpl, const int shift,
		uint64_t query_int, int c, const int ori, int *weight,
		const int max_len) {
	uint64_t query = 0;
	int n_reads = 0;
	bwa_seq_t *branch_seq = NULL, *main_seq = NULL;
	query = shift_bit(query_int, c, hm->o->k, ori);
	branch_seq = try_short_tpl_byte_ext(hm, query, c, ori, max_len / 2);
	p_ctg_seq("BRANCH", branch_seq);

	if (ori) {
		main_seq = cut_edge_tail(main_tpl, shift, max_len / 2, ori);
		p_ctg_seq("MAIN", main_seq);
		n_reads
				= find_junction_reads(hm, branch_seq, main_seq, max_len, weight);
	} else {
		main_seq = cut_edge_tail(main_tpl, shift + hm->o->k, max_len / 2, ori);
		p_ctg_seq("MAIN", main_seq);
		n_reads
				= find_junction_reads(hm, main_seq, branch_seq, max_len, weight);
	}
	bwa_free_read_seq(1, main_seq);
	bwa_free_read_seq(1, branch_seq);
	return n_reads;
}

int val_short_tpl(hash_map *hm, uint64_t query_int, int max_c, int second_c,
		const int ori) {
	uint64_t main_query = 0, main_branch = 0;
	uint64_t second_query = 0, second_branch = 0;
	main_query = shift_bit(query_int, max_c, hm->o->k, ori);
	main_branch = try_short_tpl_ext(hm, main_query, ori, SHORT_BRANCH_SHIFT);
	second_query = shift_bit(query_int, second_c, hm->o->k, ori);
	second_branch
			= try_short_tpl_ext(hm, second_query, ori, SHORT_BRANCH_SHIFT);
	if (main_branch == second_branch)
		return 0;
	// If the second branch cannot extend to SHORT_BRANCH_LEN, just prune it
	if (second_branch == LOWER_ONES_32)
		return 0;
	// If the main branch cannot extend to SHORT_BRANCH_LEN, try branching as well
	return 1;
}

int right_connect(edge *branch, hash_map *hm, tpl_hash *all_tpls,
		uint64_t query_int, const int ori) {
	int *counters = NULL, locus = 0, i = 0, connected = 0;
	uint64_t eg_id = 0, value = 0;
	edge *right_tpl = NULL;
	counters = count_next_kmers(hm, query_int, 0, ori);
	for (i = 0; i < 4; i++) {
		if (counters[i] == 0)
			continue;
		query_int = shift_bit(query_int, i, hm->o->k, ori);
		read_tpl_using_kmer(query_int, hm, &eg_id, &locus, &value);
		tpl_hash::iterator it = all_tpls->find(query_int);
		if (it != all_tpls->end()) {
			right_tpl = it->second;
			show_debug_msg(__func__,
					"Right connect [%d, %d] to [%d, %d] at %d. \n", branch->id,
					branch->len, right_tpl->id, right_tpl->len, locus);
			// This locus is not correct if connect to itself.
			add_a_junction(right_tpl, branch, query_int, locus, ori,
					counters[i]);
			connected = 1;
		}
	}
	free(counters);
	return connected;
}

/**
 * For an existing edge, try to assemble the branches
 */
void kmer_ext_branch(edge *eg, hash_map *hm, tpl_hash *all_tpls, const int ori) {
	int *counters = NULL, weight = 0, c = 0, i = 0, j = 0;
	int kmer_len = hm->o->k, branch_is_valid = 0;
	uint64_t query_int = 0, branch_query = 0;
	edge *branch_eg = NULL;
	if (eg->len <= kmer_len)
		return;
	show_debug_msg(__func__, "Branching [%d, %d]...\n", eg->id, eg->len);
	for (i = 0; i < eg->len - kmer_len; i++) {
		if (ori) {
			query_int = get_kmer_int(eg->contig->seq, eg->len - kmer_len - i,
					1, kmer_len);
		} else
			query_int = get_kmer_int(eg->contig->seq, i, 1, kmer_len);
		counters = count_next_kmers(hm, query_int, 1, ori);

		bwa_seq_t *debug = get_kmer_seq(query_int, 25);
		p_query(__func__, debug);
		bwa_free_read_seq(1, debug);
		show_debug_msg(__func__, "%d pos counters ori %d: [%d, %d, %d, %d]\n",
				i, ori, counters[0], counters[1], counters[2], counters[3]);

		for (j = 0; j < 4; j++) {
			c = counters[j];
			if (c <= 0)
				continue;
			weight = 0;
			branch_query = shift_bit(query_int, j, kmer_len, ori);
			branch_is_valid = val_branching(hm, eg, i, query_int, j, ori,
					&weight, (hm->o->read_len - SHORT_BRANCH_SHIFT) * 2);
			if (branch_is_valid) {
				branch_eg = blank_edge(branch_query, kmer_len, 1, ori);
				set_tail(eg, i, hm->o->read_len - SHORT_BRANCH_SHIFT, ori);
				kmer_ext_edge(branch_eg, branch_query, hm, all_tpls, ori);
				if (branch_eg->len >= MIN_BRANCH_LEN) {
					g_mutex_lock(kmer_id_mutex);
					if (ori)
						add_a_junction(eg, branch_eg, query_int, i, ori, weight);
					else
						add_a_junction(eg, branch_eg, query_int, i + kmer_len,
								ori, weight);
					all_tpls->insert(make_pair<uint64_t, edge*> (branch_eg->id,
							branch_eg));
					g_mutex_unlock(kmer_id_mutex);
					// Try to extend branches of current branch
					mark_tpl_kmers_used(eg, hm, kmer_len);
					upd_tpl_jun_locus(eg, branching_events, kmer_len);
					kmer_ext_branch(branch_eg, hm, all_tpls, ori);
				} else {
					mark_tpl_kmers_fresh(eg, hm, kmer_len);
					destroy_eg(branch_eg);
				}
			}
		}
		free(counters);
		counters = NULL;
	}
	free(counters);
}

void kmer_ext_edge(edge *eg, uint64_t query_int, hash_map *hm,
		tpl_hash *all_tpls, const int ori) {
	int max_c = 0, second_c = 0, *counters = NULL, weight = 0;
	int main_is_valid = 0, branch_is_valid = 0, short_val_rs = 0;
	bwa_seq_t *debug = NULL;

	show_debug_msg(__func__, "------ Started extending edge %d to ori %d...\n",
			eg->id, ori);
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	while (1) {
		weight = 0;
		counters = count_next_kmers(hm, query_int, 1, ori);
		max_c = get_max_index(counters);
		second_c = get_second_freq_char(counters, max_c);

		if (ori)
			seq_reverse(eg->len, eg->contig->seq, 0);
		debug = get_kmer_seq(query_int, 25);
		p_query(__func__, debug);
		bwa_free_read_seq(1, debug);
		show_debug_msg(__func__, "Next chars: [%d, %d, %d, %d] \n",
				counters[0], counters[1], counters[2], counters[3]);
		show_debug_msg(__func__,
				"Ori %d, Edge %d, length %d, Next char: %d \n", ori, eg->id,
				eg->len, max_c);
		p_ctg_seq("Contig", eg->contig);
		if (ori)
			seq_reverse(eg->len, eg->contig->seq, 0);

		if (max_c == -1) {
			if (right_connect(eg, hm, all_tpls, query_int, ori)) {
				show_debug_msg(__func__, "[%d, %d] Right connected. \n",
						eg->id, eg->len);
			} else
				show_debug_msg(__func__, "[%d, %d] No hits, stop here. \n",
						eg->id, eg->len);
			break;
		}
		// The main branch may be not valid, check the second frequent char, try it
		if (second_c != -1 && counters[second_c] > 1) {
			//show_debug_msg(__func__, "[%d, %d] Validating short edges... \n",
			//		eg->id, eg->len);
			short_val_rs = val_short_tpl(hm, query_int, max_c, second_c, ori);
			if (short_val_rs) {
				main_is_valid = val_junction_reads(hm, eg, query_int, max_c,
						ori, &weight, (hm->o->read_len - SHORT_BRANCH_SHIFT)
								* 2);
				// If the main branch is not valid, set the next char to second frequent char.
				if (!main_is_valid) {
					branch_is_valid = val_junction_reads(hm, eg, query_int,
							second_c, ori, &weight, (hm->o->read_len
									- SHORT_BRANCH_SHIFT) * 2);
					if (branch_is_valid)
						max_c = second_c;
					else {
						free(counters);
						break;
					}
				}
			}
		}
		free(counters);
		ext_con(eg->contig, max_c, 0);
		eg->len = eg->contig->len;
		// In case the template has repeats, so mark kmers used TEMPERARIALY. The locus is incorrect.
		mark_kmer_used(query_int, hm, eg->id, eg->len);
		query_int = shift_bit(query_int, max_c, hm->o->k, ori);
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, Edge %d, length %d \n", ori,
					eg->id, eg->len);
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
}

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	edge *eg = NULL;
	int round_1_len = 0, round_2_len = 0;
	uint64_t kmer_int = 0, rev_kmer_int = 0;
	map_opt *opt = NULL;
	kmer_counter *counter = NULL;
	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	tpl_hash *all_tpls = params->all_tpls;
	opt = params->hm->o;

	counter = (kmer_counter*) data;
	kmer_int = counter->kmer;
	rev_kmer_int = rev_comp_kmer(kmer_int, opt->k);
	if (kmer_int < 64 || rev_kmer_int < 64 || get_kmer_count(kmer_int,
			params->hm, 1) <= 1) {
		return NULL;
	}
	eg = blank_edge(kmer_int, opt->k, opt->k, 0);

	show_debug_msg(__func__, "============= %" ID64 ": %d ============ \n",
			kmer_int, counter->count);
	kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);

	kmer_int = get_kmer_int(eg->contig->seq, 0, 1, opt->k);
	kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 1);
	show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		kmer_int = get_kmer_int(eg->contig->seq, eg->len - opt->k, 1, opt->k);
		round_2_len = eg->len;
		kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 0);
		show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			kmer_int = get_kmer_int(eg->contig->seq, 0, 1, opt->k);
			kmer_ext_edge(eg, kmer_int, params->hm, all_tpls, 1);
			show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id,
					eg->len);
		}
	}

	mark_tpl_kmers_used(eg, params->hm, opt->k);
	upd_tpl_jun_locus(eg, branching_events, opt->k);
	kmer_ext_branch(eg, params->hm, all_tpls, 0);
	kmer_ext_branch(eg, params->hm, all_tpls, 1);

	eg->start_kmer_int = *((uint64_t*) data);
	if (eg->len <= opt->k) {
		destroy_eg(eg);
	} else {
		g_mutex_lock(kmer_id_mutex);
		all_tpls->insert(make_pair<uint64_t, edge*> (eg->id, eg));
		g_mutex_unlock(kmer_id_mutex);
	}
	return NULL;
}

void kmer_threads(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0, *freq = 0, count = 0;
	hash_map *hm = params->hm;
	mer_hash *hash = hm->hash;
	kmer_counter *counter = NULL;
	GPtrArray *start_kmers = g_ptr_array_sized_new(BUFSIZ);

	for (mer_hash::iterator it = hash->begin(); it != hash->end(); ++it) {
		freq = it->second;
		count = freq[0];
		if (count > 1 && it->first > 64 && rev_comp_kmer(it->first, hm->o->k)
				> 64) {
			counter = (kmer_counter*) malloc(sizeof(kmer_counter));
			counter->kmer = it->first;
			counter->count = count;
			g_ptr_array_add(start_kmers, counter);
		}
	}
	show_msg(__func__, "Sorting %d initial kmers... \n", start_kmers->len);
	g_ptr_array_sort(start_kmers, (GCompareFunc) cmp_kmers_by_count);
	show_msg(__func__, "Extending by kmers...\n");
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);
	for (i = 0; i < start_kmers->len; i++) {
		counter = (kmer_counter*) g_ptr_array_index(start_kmers, i);
		g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		//if (i > 11000)
		//	break;
	}
	for (i = 0; i < start_kmers->len; i++) {
		counter = (kmer_counter*) g_ptr_array_index(start_kmers, i);
		free(counter);
	}
	g_ptr_array_free(start_kmers, TRUE);
	g_thread_pool_free(thread_pool, 0, 1);
}

void store_junctions(char *name, GPtrArray *branching_events) {
	junction *jun = NULL;
	uint64_t i = 0;
	char entry[BUFSIZE];
	FILE *f = xopen(name, "w");
	sprintf(entry, "Main\tBranch\tLocus\tWeight\tDirection\n");
	fputs(entry, f);
	for (i = 0; i < branching_events->len; i++) {
		jun = (junction*) g_ptr_array_index(branching_events, i);
		sprintf(entry, "[%d, %d]\t[%d, %d]\t%d\t%d\t%d\n", jun->main_tpl->id,
				jun->main_tpl->len, jun->branch_tpl->id, jun->branch_tpl->len,
				jun->locus, jun->weight, jun->ori);
		fputs(entry, f);
	}
	fclose(f);
}

void test_kmer_ext(kmer_t_meta *params) {
	uint64_t kmer_int = 896094812036308;
	int round_1_len = 0, round_2_len = 0;
	edge *eg = new_eg();
	FILE *contigs = NULL;
	GPtrArray *all_edges = NULL;
	map_opt *opt = params->hm->o;
	tpl_hash *all_tpls = params->all_tpls;

	eg = blank_edge(kmer_int, 25, 25, 0);
	kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 0);
	round_1_len = eg->len;
	show_debug_msg(__func__, "Extending to the right: [%d, %d]... \n", eg->id,
			eg->len);
	kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 1);
	show_debug_msg(__func__, "Done Extending to the right: [%d, %d]... \n",
			eg->id, eg->len);

	//	if (eg->len - round_1_len > 2) {
	//		kmer_int = get_kmer_int(eg->contig->seq, eg->len - opt->k, 1, opt->k);
	//		round_2_len = eg->len;
	//		kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 0);
	//		// Add the new length to existing junctions
	//		update_junctions(branching_events, eg, eg->len - round_2_len);
	//		show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id, eg->len);
	//		if (eg->len - round_2_len > 2) {
	//			kmer_int = get_kmer_int(eg->contig->seq, 0, 1, opt->k);
	//			kmer_ext_edge(eg, kmer_int, params->hm, params->all_tpls, 1);
	//			show_debug_msg(__func__, "Edge %d with length: %d \n", eg->id,
	//					eg->len);
	//		}
	//		update_junctions(branching_events, eg, 0);
	//	}

	kmer_ext_branch(eg, params->hm, all_tpls, 0);
	kmer_ext_branch(eg, params->hm, all_tpls, 1);
	all_tpls->insert(make_pair<uint64_t, edge*> (eg->id, eg));

	all_edges = hash_to_array(all_tpls);
	contigs = xopen(get_output_file("single.fa", kmer_out), "w");
	save_edges(all_edges, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);
	store_junctions(get_output_file("single.junctions", kmer_out),
			branching_events);
}

void ext_by_kmers_core(char *lib_file, const char *solid_file) {
	hash_map *hm = NULL;
	mer_hash map;
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *kmer_edges = NULL;

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	hm = load_hash_map(lib_file, 0, map);
	hm->hash = &map;

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers hash map: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	branching_events = g_ptr_array_sized_new(BUFSIZ);
	params->hm = hm;
	params->all_tpls = &all_tpls;

	test_kmer_ext(params);
	exit(1);
	kmer_threads(params);
	//pick_unused_kmers(params);

	contigs = xopen(get_output_file("paired.fa", kmer_out), "w");
	kmer_edges = hash_to_array(&all_tpls);
	save_edges(kmer_edges, contigs, 0, 0, 0);
	fflush(contigs);
	fclose(contigs);

	store_junctions(get_output_file("paired.junctions", kmer_out),
			branching_events);
	//	reset_edge_ids(all_kmer_edges);
	//	merge_ol_edges(all_kmer_edges, ins_size, sd_ins_size, hm->seqs,
	//			kmer_n_threads);
	//	free(params);
	//	destroy_hm(hm);
	//
	//	contigs = xopen(get_output_file("merged.fa", kmer_out), "w");
	//	save_edges(all_kmer_edges, contigs, 0, 0, 100);
	//	fflush(contigs);
	//	fclose(contigs);
	//	contigs = xopen(get_output_file("merged_all.fa", kmer_out), "w");
	//	save_edges(all_kmer_edges, contigs, 0, 0, 0);
	//	fflush(contigs);
	//	fclose(contigs);
}

int pe_kmer(int argc, char *argv[]) {
	int c = 0;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);
	while ((c = getopt(argc, argv, "k:m:s:o:t:")) >= 0) {
		switch (c) {
		case 'k':
			kmer_len = atoi(optarg);
			break;
		case 'o':
			kmer_out = optarg;
			break;
		case 'm':
			ins_size = atoi(optarg);
			break;
		case 's':
			sd_ins_size = atoi(optarg);
			break;
		case 't':
			kmer_n_threads = atoi(optarg);
			break;
		}
	}
	if (optind + 2 > argc) {
		show_msg(__func__, "Parameters error! \n");
		return 1;
	}
	if (!g_thread_supported())
		g_thread_init(NULL);
	kmer_id_mutex = g_mutex_new();

	ext_by_kmers_core(argv[optind], argv[optind + 1]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
