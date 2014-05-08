/*
 * ass.cpp
 *
 *  Created on: 09-Apr-2013
 *      Author: carl
 */

#include <unordered_map>
#include <stdint.h>
#include <glib.h>
#include <cstdio>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bwtaln.h"
#include "tpl.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "ass.hpp"
#include "junction.hpp"
#include "graph.hpp"
#include "pool.hpp"
#include "merge.hpp"
#include "hash.hpp"
#include "psl.h"

using namespace std;

int TESTING = 1563032;
int TESTING2 = 514485;
int MAX_FRESH_TRIAL = 2;
int DETAIL_ID = -1;

int test_suffix = 0;
int kmer_ctg_id = 1;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
int stage = 1;
int fresh_trial = 0;
uint32_t n_used_reads = 0;
char *kmer_out = NULL;
GPtrArray *hang_reads = NULL;

uint32_t *tpl_kmers = NULL;
uint32_t tpl_kmers_hash_size = 1;

bwa_seq_t *TEST = NULL;
bwa_seq_t *CLUSTER_START_READ = NULL;

GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

void reset_tpl_hash() {
	if(tpl_kmers) free(tpl_kmers);
	tpl_kmers = (uint32_t*) calloc(tpl_kmers_hash_size, sizeof(uint32_t));
}

/**
 * Initialize a template
 */
tpl *blank_tpl(bwa_seq_t *start_read, int len, int ori, char *step) {
	tpl *t = new_tpl();
	g_mutex_lock(kmer_id_mutex);
	t->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	if (start_read->rev_com)
		switch_fr(start_read);
	if (len == 0) {
		t->ctg = blank_seq(start_read->len);
	} else {
		if (len < 0 || len > start_read->len) {
			t->ctg = new_seq(start_read, start_read->len, 0);
		} else {
			t->ctg = ori ? new_seq(start_read, len, 0) : new_seq(start_read,
					len, start_read->len - len);
		}
	}
	if (start_read->rev_com)
		switch_fr(start_read);
	t->start_read = start_read;//CLUSTER_START_READ;
	t->len = t->ctg->len;
	t->ctg->rev_com = 0;
	t->step = (char*) malloc(sizeof(char) * 32);
	sprintf(t->step, "%s", step);
	return t;
}

GPtrArray *hash_to_array(tpl_hash *all_tpls) {
	tpl *t = NULL;
	GPtrArray *tpls = g_ptr_array_sized_new(all_tpls->size());
	show_msg(__func__, "Putting hashed templates to array ... \n");
	for (tpl_hash::iterator m = all_tpls->begin(); m != all_tpls->end(); ++m) {
		t = (tpl*) m->second;
		//p_tpl(t);
		if (t->alive) {
			//show_debug_msg(__func__, "LAST READ [%d, %d] \n", t->id, t->len);
			//p_query("LAST READ", t->last_read);
			//p_tpl_reads(t);
			//p_junctions(t->m_juncs);
			//p_junctions(t->b_juncs);
			g_ptr_array_add(tpls, t);
		} else {
			destory_tpl_junctions(t);
			destroy_tpl(t, TRIED);
		}
	}
	g_ptr_array_sort(tpls, (GCompareFunc) cmp_tpl_by_id);
	return tpls;
}

/**
 * Use read-length tail to search,
 * 	find those templates could be connected to current branch
 */
void p_tpl_junctions(tpl_hash *all_tpls, int id) {
	tpl_hash::iterator it = all_tpls->find(id);
	if (it == all_tpls->end()) {
		show_debug_msg(__func__, "No template id %d \n", id);
		return;
	}
	tpl *t = (tpl*) it->second;
	printf("\n");
	p_tpl(t);
	show_debug_msg(__func__, "Template [%d, %d]: as main branches: \n", t->id,
			t->len);
	p_junctions(t->m_juncs);
	show_debug_msg(__func__, "Template [%d, %d]: as branch branches: \n",
			t->id, t->len);
	p_junctions(t->b_juncs);
	show_debug_msg(__func__, "--- End of printing junctions. --- \n\n");
}

/**
 * Remove a dead template from the global template hash
 */
void rm_global_tpl(tpl_hash *all_tpls, tpl *t, int status) {
	if (!t->alive) {
		show_debug_msg(__func__, "Template [%d, %d] is destroyed.\n", t->id,
				t->len);
		g_mutex_lock(kmer_id_mutex);
		all_tpls->erase(t->id);
		g_mutex_unlock(kmer_id_mutex);
		destory_tpl_junctions(t);
		destroy_tpl(t, status);
	}
}

tpl *add_global_tpl(tpl_hash *all_tpls, bwa_seq_t *branch_read, char *step,
		int len, int ori) {
	tpl *branch = blank_tpl(branch_read, len, ori, step);
	g_mutex_lock(kmer_id_mutex);
	all_tpls->insert(make_pair<int, tpl*> ((int) branch->id, (tpl*) branch));
	g_mutex_unlock(kmer_id_mutex);
	return branch;
}

void p_test_read() {
	if (TESTING)
	p_query("TEST", TEST);
	//	if (TEST->status == FRESH && TEST->pos != IMPOSSIBLE_NEGATIVE) {
	//		show_debug_msg(__func__, "ID: %d \n", kmer_ctg_id);
	//		exit(1);
	//	}
}

/**
 * Mark the reads on a template as TRIED temporarily
 */
void mark_as_hang_tmp(tpl *t) {
	bwa_seq_t *r = 0;
	int i = 0;
	if (t->reads) {
		for (i = 0; i < t->reads->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
			reset_to_hang(r);
			g_ptr_array_add(hang_reads, r);
		}
	}
	if (t->tried) {
		for (i = 0; i < t->tried->len; i++) {
			r = (bwa_seq_t*) g_ptr_array_index(t->tried, i);
			reset_to_hang(r);
			g_ptr_array_add(hang_reads, r);
		}
	}
}

void copy_to_hang(GPtrArray *reads) {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(reads, i);
		g_ptr_array_add(hang_reads, r);
		reset_to_hang(r);
	}
}

void set_hang_status() {
	int i = 0;
	bwa_seq_t *r = NULL;
	for (i = 0; i < hang_reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hang_reads, i);
		reset_to_hang(r);
	}
}

/**
 * Unfrozen all TRIED reads
 */
void unfrozen_hang_reads() {
	bwa_seq_t *r = NULL;
	while (hang_reads->len > 0) {
		r = (bwa_seq_t*) g_ptr_array_index(hang_reads, 0);
		//if (r->status == HANG)
		reset_to_fresh(r);
		g_ptr_array_remove_index_fast(hang_reads, 0);
	}
}

/**
 * Extend a template until no next kmer
 */
int kmer_ext_tpl(hash_table *ht, tpl_hash *all_tpls, pool *p, tpl *from,
		tpl *t, bwa_seq_t *query, const int ori) {
	int max_c = -1;
	int to_connect = 0;
	int ext_len = 0, no_read_len = 0;
	GPtrArray *near_tpls = g_ptr_array_sized_new(4);
	bwa_seq_t *tail = new_seq(query, query->len, 0), *adj_tail = NULL;
	bwa_seq_t *last_read = NULL;

	show_debug_msg(__func__, "------ Extending tpl %d to ori %d ... ------\n", t->id, ori);
	p_query(__func__, tail);
	near_tpls = nearby_tpls(t, 1);
	if (from) g_ptr_array_add(near_tpls, from);
	p_test_read();
	while (1) {
		// If the query is bad, stop
		if (is_bad_query(tail)) {
			show_debug_msg(__func__, "Repetitive tail, stop.\n");
			p_query(__func__, tail);
			mark_pool_reads_tried(p, t);
			to_connect = 0;
			break;
		}

		if (p->reads->len > 0)
			last_read = (bwa_seq_t*) g_ptr_array_index(p->reads, p->reads->len - 1);

		if (t->id == DETAIL_ID || DETAIL_ID == 0) {
			p_test_read();
			p_query("QUERY", tail);
			p_ctg_seq("TEMPLATE", t->ctg);
			p_pool("CURRENT POOL", p, NULL);
		}

		max_c = get_next_char(ht, p, t, ori);

		if (t->id == DETAIL_ID || DETAIL_ID == 0)
			show_debug_msg(__func__,
					"Ori: %d, Template [%d, %d], Next char: %c \n", ori, t->id,
					t->len, "ZACGTN"[max_c + 1]);

		// If cannot extend or the reads in the pool is too few, try to add mates into the pool
		if (max_c == -1) {
			//p_ctg_seq("TEMPLATE", t->ctg);
			//show_debug_msg(__func__, "Looking for mates on [%d, %d] ...\n", t->id, t->len);
			//p_test_read();
			//p_tpl_reads(t);
			find_ol_mates(ht, p, near_tpls, t, tail, N_MISMATCHES, ori);
			//p_pool("MATE_POOL", p, NULL);
			max_c = get_next_char(ht, p, t, ori);
			if (max_c == -1) {
				if (adj_tail) {
					bwa_free_read_seq(1, tail);
					tail = adj_tail;
					max_c = get_next_char(ht, p, t, ori);
				}
				if (max_c == -1) {
					show_debug_msg(__func__, "No hits, stop ori %d: [%d, %d] \n", ori, t->id, t->len);
					if (last_read) t->last_read = last_read;
					to_connect = 1; break;
				}
			} else {
				p_query("TAIL", tail);
				show_debug_msg(__func__, "Added mates: ori %d \n", ori);
				p_pool("MATE_POOL", p, NULL);
			}
		}

		ext_con(t->ctg, max_c, ori);
		t->len = t->ctg->len;
		ext_que(tail, max_c, ori);
		ext_len++;
		// If the overlapped region between t and r has too many mismatches, remove from pool
		rm_half_clip_reads(p, t, max_c, N_MISMATCHES, ori);
		if (forward(p, t, ori)) no_read_len = 0;
		else no_read_len++;
		next_pool(ht, p, near_tpls, t, tail, N_MISMATCHES, ori);

		if (t->len % 100 == 0)
			show_debug_msg(__func__, "Ori %d, tpl %d, length %d \n", ori, t->id, t->len);
		if (t->len > 100000) {
			p_tpl(t);
			err_fatal(__func__,
					"[ERROR] Too long contig [%d, %d]. Please report to caishaojiang@gmail.com.\n",
					t->id, t->len);
		}
	}
	set_rev_com(t->ctg);
	g_ptr_array_free(near_tpls, TRUE);
	bwa_free_read_seq(1, tail);
	p_test_read();
	return to_connect;
}

/**
 * This is a unit to extend a template to left/right.
 * @to_connect: if 1, try to connect to existing templates after extension
 * @from: if doing jumping, we jump from some template 'from', whose reads are used
 * 			for mate pool determination
 * @init_p: initial pool; if NULL, initialise the pool
 * @init_q: initial query; if NULL, get template tail as query
 */
int ext_unit(hash_table *ht, tpl_hash *all_tpls, pool *init_p, tpl *from,
		tpl *t, bwa_seq_t *init_q, int to_connect, int ori) {
	pool *p = init_p;
	int to_con = 0, connected = 0;
	bwa_seq_t *query = NULL;

	if (!t || !t->alive) return -1;

	if (!init_p) {
		p = new_pool();
		init_pool(ht, p, t, kmer_len, N_MISMATCHES, ori);
		query = get_tail(t, kmer_len, ori);
	} else {
		query = init_q;
	}
	if (p->reads->len == 0 && t->len <= ht->o->read_len) {
		bwa_free_read_seq(1, query);
		destroy_pool(p);
		return -1;
	}
	// Correct the initial bases if the starting read is not good
	if (t->len == ht->o->read_len)
		correct_init_tpl_base(p, t, ori);
	// If extending to left, reverse the locus first
	if (ori) reverse_locus(t);
	to_con = kmer_ext_tpl(ht, all_tpls, p, from, t, query, ori);
	bwa_free_read_seq(1, query);
	destroy_pool(p);

	if (!t->alive) return -1;
	set_rev_com(t->ctg);
	if (ori) upd_locus_on_tpl(t, 0, 0);
	return connected;
}

/**
 * Extend a read to a template;
 * No branching and validation
 */
tpl *ext_a_read(hash_table *ht, tpl_hash *all_tpls, tpl *from, bwa_seq_t *read, index64 count) {
	tpl *t = NULL;
	int iter = 0, after_unit = -1, pre_len = 0;
	if (is_biased_q(read) || has_n(read, 1) || is_repetitive_q(read) || read->status != FRESH) {
		return NULL;
	}
	printf("\n");
	show_debug_msg(__func__, "============= FRESH %s: %" ID64 " ============ \n", read->name, count);

	fresh_trial++;
	p_query(__func__, read);
	t = add_global_tpl(all_tpls, read, "FRESH", read->len, 0);
	t->kmer_freq = count;
	if (count > 1) mark_init_reads_used(ht, t, read, N_MISMATCHES);
	else add2tpl(t, read, 0);
	// Right->left->right->left...until not extendable
	while (iter++ <= 4 && t->len > pre_len) {
		// Extend to the left first
		after_unit = ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 1);
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
		if (iter > 1 && after_unit == -1) break;
		pre_len = t->len;
		after_unit = ext_unit(ht, all_tpls, NULL, from, t, NULL, 0, 0);
		show_debug_msg(__func__, "tpl %d with length: %d \n", t->id, t->len);
		if (after_unit == -1) break;
	}
	// If the read cannot be even extend one base, reset the read to fresh
	if (t->len == read->len) {
		reset_to_fresh(read);
	}
	return t;
}

void load_anchors(tpl_hash *all_tpls, hash_table *tpl_ht, int8_t *candiates, GPtrArray *anchors, tpl *t, int rev_com, int s, int e) {
	int i = 0, j = 0, locus = 0, start = 0, end = 0;
	hash_key key = 0ULL;
	index64 tpl_id = 0;
	hash_value value = 0ULL;
	anchor *a = NULL;
	tpl *sharing = NULL;
	bwa_seq_t *kmer = NULL;
	//p_hash_table(tpl_ht);
	//show_debug_msg(__func__, "Templates sharing 11-mer with template [%d, %d]; rev_com: %d \n", t->id, t->len, rev_com);
	for (i = max(0, s); i <= min(t->len, e) - tpl_ht->o->k; i++) {
		key = get_hash_key(t->ctg->seq, i, tpl_ht->o->interleaving, tpl_ht->o->k);
		if (rev_com) key = get_hash_key(t->ctg->rseq, i, tpl_ht->o->interleaving, tpl_ht->o->k);
		start = tpl_ht->k_mers_occ_acc[key];
		end = (key >= tpl_ht->o->n_k_mers) ? tpl_ht->k_mers_occ_acc[tpl_ht->o->n_k_mers - 1]
				: tpl_ht->k_mers_occ_acc[key + 1];
		if (end - start >= 1) {
			//kmer = get_key_seq(key, tpl_ht->o->k);
			//kmer->contig_locus = i; kmer->contig_id = t->id;
			//p_query("KMER", kmer);
			//bwa_free_read_seq(1, kmer);
			for (j = start; j < end; j++) {
				value = tpl_ht->pos[j];
				read_hash_value(&tpl_id, &locus, value);
				if (tpl_id == t->id || candiates[tpl_id] != 1) continue;
				//show_debug_msg(__func__, "HASH = %"ID64"; At %d: template %d at %d \n", value, i, tpl_id, locus);
				tpl_hash::iterator it = all_tpls->find(tpl_id);
				if (it == all_tpls->end()) continue;
				sharing = (tpl*) it->second;
				if (sharing->alive) {
					a = (anchor*) malloc(sizeof(anchor));
					a->b = sharing;
					a->size = tpl_ht->o->k;
					a->from = i;
					a->locus = locus;
					a->hash = value;
					a->rev_com = rev_com;
					g_ptr_array_add(anchors, a);
					//p_anchor("HIT", a);
				}
			}
		}
	}
}

void tpls_sharing_kmers(hash_table *ht, tpl_hash *all_tpls, hash_table *tpl_ht, GPtrArray *anchors,
		tpl *t, int s, int e, int rev_com) {
	int i = 0;
	bwa_seq_t *kmer = NULL, *r = NULL, *m = NULL;
	anchor *a = NULL, *pre = NULL;
	if (!t->alive || t->reads->len <= 0 || t->pair_pc >= 1.0) return;
	// It is a list, acting as a quick hash table, indicating the templates that have reads pairing with 't'
	// For two templates, only if they have pairs, they could be probably connected.
	int8_t *candiates = (int8_t*) calloc(kmer_ctg_id + 4, sizeof(int8_t));
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, ht->seqs);
		if (m->status == USED && m->contig_id != t->id) {
			//p_query("READ", r);
			//p_query("MATE", m);
			candiates[m->contig_id] = 1;
		}
	}
	load_anchors(all_tpls, tpl_ht, candiates, anchors, t, rev_com, s, e);
	free(candiates);
	compact_anchors(anchors);
	return;
}

void branching(hash_table *ht, tpl_hash *all_tpls, GPtrArray *tpls) {
	hash_table *tpl_ht = hash_tpls(tpls, ht->o->k, 1);
	GPtrArray *anchors = NULL;
	int i = 0, j = 0, n_both_connected = 0, n_one_connected = 0;
	int both_connected = 0, one_connected = 0;
	tpl *t = NULL;
	anchor *a = NULL;
	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		//printf("\n------------------------------------------------------- \n");
		//p_tpl(t);
		if (!t->alive || t->pair_pc >= 1.0) continue;
		anchors = g_ptr_array_sized_new(4);
		tpls_sharing_kmers(ht, all_tpls, tpl_ht, anchors, t, 0, t->len, 0);
		both_connected = connect_both_ends(ht, anchors, t);
		n_both_connected += both_connected;
		for (j = 0; j < anchors->len; j++) {
			a = (anchor*) g_ptr_array_index(anchors, j); free(a);
		}
		g_ptr_array_free(anchors, TRUE);

		// If the main and branch templates are with different orientation
		if (!both_connected) {
			anchors = g_ptr_array_sized_new(4);
			tpls_sharing_kmers(ht, all_tpls, tpl_ht, anchors, t, 0, t->len, 1);
			rev_com_tpl(t);
			both_connected = connect_both_ends(ht, anchors, t);
			if (both_connected) n_both_connected++;
			else rev_com_tpl(t);
			for (j = 0; j < anchors->len; j++) {
				a = (anchor*) g_ptr_array_index(anchors, j); free(a);
			}
			g_ptr_array_free(anchors, TRUE);
		}

		anchors = NULL;
	}
	show_msg(__func__, "%d templates are both connected. \n", n_both_connected);
	show_msg(__func__, "Perform branching ... \n");
	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		if (t->id == 1) rev_com_tpl(t);
		if (!t->alive || t->pair_pc >= 1.0) continue;
		anchors = g_ptr_array_sized_new(4);
		tpls_sharing_kmers(ht, all_tpls, tpl_ht, anchors, t, 0, t->len, 0);
		one_connected = connect_one_end(ht, anchors, t);
		n_one_connected += one_connected;
		for (j = 0; j < anchors->len; j++) {
			a = (anchor*) g_ptr_array_index(anchors, j); free(a);
		}
		g_ptr_array_free(anchors, TRUE);
	}
	show_msg(__func__, "%d branches are created. \n", n_one_connected);
	destory_tpl_ht(tpl_ht);
}

int connect_paired_tpls(kmer_t_meta *params, GPtrArray *tpls) {
	int connected = 0, both_connected = 0, n_merged = 1, n_both_connected = 0, i = 0, j = 0, iter = 1;
	tpl *t = NULL, *b = NULL;
	anchor *a = NULL;
	hash_table *ht = params->ht;
	tpl_hash *all_tpls = params->all_tpls;
	GPtrArray *anchors = NULL;
	hash_table *tpl_ht = NULL;
	int right_tpl_id = 0;
	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		t->visited = 0;
		mv_paired_reads_to_tried(ht, t);
	}
	while(iter < 16 && n_merged) {
		n_merged = 0; n_both_connected = 0;
		tpl_ht = hash_tpls(tpls, ht->o->k, 1);
		//p_hash_table(tpl_ht);
		show_msg(__func__, "Connecting %d templates, iteration %d ... \n", tpls->len, iter++);
		for (i = 0; i < tpls->len; i++) {
			t = (tpl*) g_ptr_array_index(tpls, i);
			//p_tpl(t);
			if (!t->alive || t->pair_pc >= 1.0) continue;
			connected = 0;

			// Both template are with the same orientation
			anchors = g_ptr_array_sized_new(4);
			tpls_sharing_kmers(ht, params->all_tpls, tpl_ht, anchors, t, 0, t->len, 0);
			for (j = 0; j < anchors->len; j++) {
				a = (anchor*) g_ptr_array_index(anchors, j);
				//p_anchor("HIT", a);
				if (!connected && a->size > 0 && a->b->alive && !a->b->visited) {
					connected = connect_at_locus_right(ht, t, a->b, a->from, a->locus);
					if (connected)  {
						t->visited = 1; n_merged++;
					}
				}
				free(a);
			}
			g_ptr_array_free(anchors, TRUE);
			if (connected) continue;

			// If the two templates are not with the same orientation
			anchors = g_ptr_array_sized_new(4);
			rev_com_tpl(t);
			tpls_sharing_kmers(ht, params->all_tpls, tpl_ht, anchors, t, 0, t->len, 0);
			for (j = 0; j < anchors->len; j++) {
				a = (anchor*) g_ptr_array_index(anchors, j);
				//p_anchor("HIT", a);
				if (!connected && a->size > 0 && a->b->alive && !a->b->visited) {
					connected = connect_at_locus_right(ht, t, a->b, a->from, a->locus);
					if (connected)  {
						t->visited = 1; n_merged++; break;
					}
				}
			}
			// Reverse complement the right template
			if (!connected) {
				rev_com_tpl(t);
				for (j = 0; j < anchors->len; j++) {
					a = (anchor*) g_ptr_array_index(anchors, j);
					//p_anchor("HIT", a);
					if (!connected && a->size > 0 && a->b->alive && !a->b->visited) {
						rev_com_tpl(a->b);
						connected = connect_at_locus_right(ht, t, a->b, t->len - a->from - a->size,
								a->b->len - a->locus - a->size);
						if (connected)  {
							t->visited = 1; n_merged++; break;
						} else rev_com_tpl(a->b);
					}
				}
			}
			for (j = 0; j < anchors->len; j++) {
				a = (anchor*) g_ptr_array_index(anchors, j);
				free(a);
			}
			g_ptr_array_free(anchors, TRUE);

			// If the are enough spanning pairs, merge them
			if (!connected) {
				right_tpl_id = right_tpl_to_merge(ht->seqs, t, PAIR_PERCENTAGE);
				show_debug_msg(__func__, "Right template of [%d, %d]: %d\n", t->id, t->len, right_tpl_id);
				if (right_tpl_id >= 0) {
					tpl_hash::iterator it = all_tpls->find(right_tpl_id);
					if (it != all_tpls->end()) {
						b = (tpl*) it->second;
						if (left_tpl_is_paired(ht->seqs, t, b, PAIR_PERCENTAGE))
							connected = merged_jumped(ht, t, b, N_MISMATCHES);
					}
				}
				if (connected) n_merged++;
				// Reverse complement the current left template
				if (!connected) {
					rev_com_tpl(t);
					right_tpl_id = right_tpl_to_merge(ht->seqs, t, PAIR_PERCENTAGE);
					show_debug_msg(__func__, "Right template of reverse [%d, %d]: %d\n", t->id, t->len, right_tpl_id);
					if (right_tpl_id >= 0) {
						tpl_hash::iterator it = all_tpls->find(right_tpl_id);
						if (it != all_tpls->end()) {
							b = (tpl*) it->second;
							if (left_tpl_is_paired(ht->seqs, t, b, PAIR_PERCENTAGE))
								connected = merged_jumped(ht, t, b, N_MISMATCHES);
						}
					}
					if (connected)  {
						t->visited = 1; n_merged++;
					} else rev_com_tpl(t);
				}
			}
		}
		for (i = 0; i < tpls->len; i++) {
			t = (tpl*) g_ptr_array_index(tpls, i);
			t->visited = 0;
		}
		show_msg(__func__, "%d templates are merged. \n", n_merged);
		if (n_merged > 0) {
			destory_tpl_ht(tpl_ht); tpl_ht = NULL;
		}
	}
	// Connect both end of a branch template to a main template
	if (!tpl_ht) destory_tpl_ht(tpl_ht);
	branching(ht, all_tpls, tpls);
	int n_pairs = 0;
	for (i = 0; i < tpls->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls, i);
		if (!t->alive) {
			rm_global_tpl(params->all_tpls, t, FRESH);
			g_ptr_array_remove_index_fast(tpls, i--);
		} else {
			mv_paired_reads_back(t);
			n_pairs = count_pairs_on_tpl(t->reads);
			t->pair_pc = (((float) n_pairs) * 2.0) / ((float) t->reads->len);
		}
	}
	return connected;
}
/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	uint64_t read_id = 0;
	occ_counter *counter = NULL;
	bwa_seq_t *read = NULL;
	tpl *t = NULL;

	kmer_t_meta *params = (kmer_t_meta*) thread_params;
	tpl_hash *all_tpls = params->all_tpls;
	hash_table *ht = params->ht;

	counter = (occ_counter*) data;
	read_id = counter->kmer;
	read = &ht->seqs[read_id];

	if (counter->count < 1)  return NULL;
	if (TESTING && fresh_trial == 0) read = &ht->seqs[TESTING];
	if (TESTING && fresh_trial == 1) read = &ht->seqs[TESTING2];

	t = ext_a_read(ht, all_tpls, NULL, read, counter->count);
	if (t) {
		if (t->alive && t->len > read->len + 2) {
			unfrozen_tried(t);
			refresh_tpl_reads(ht, t, 0, t->len, N_MISMATCHES);
			correct_tpl_base(ht->seqs, t, ht->o->read_len, 0, t->len);
			//p_tpl(t);
		} else rm_global_tpl(all_tpls, t, TRIED);
	}
	return NULL;
}

GPtrArray *sort_tpls_by_pair_pc(kmer_t_meta *params) {
	int i = 0; tpl *t = NULL;
	hash_table *ht = params->ht;
	GPtrArray *tpls_by_pair_pc = hash_to_array(params->all_tpls);
	g_ptr_array_sort(tpls_by_pair_pc, (GCompareFunc) cmp_tpl_by_rev_pair_pc);
	for (i = 0; i < tpls_by_pair_pc->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls_by_pair_pc, i);
		t->visited = 0;
		if (!t->alive) {
			t->alive = 0;
			rm_global_tpl(params->all_tpls, t, FRESH);
			g_ptr_array_remove_index_fast(tpls_by_pair_pc, i--);
		}
	}
	reset_seqs(ht->seqs, ht->n_seqs);
	return tpls_by_pair_pc;
}

void connect_tpls(kmer_t_meta *params) {
	GPtrArray *tpls_by_pair_pc = sort_tpls_by_pair_pc(params);
	connect_paired_tpls(params, tpls_by_pair_pc);
	g_ptr_array_free(tpls_by_pair_pc, TRUE);
}

int connect_missing_ends(hash_table *ht, tpl_hash *all_tpls, tpl *t) {
	GPtrArray *reads = g_ptr_array_sized_new(4);
	GPtrArray *counters = NULL;
	occ_counter *c = NULL;
	int i = 0, right_tpl_id = 0, connected = 0;
	tpl *b = NULL;
	bwa_seq_t *r = NULL, *m = NULL;
	for (i = 0; i < t->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(t->reads, i);
		m = get_mate(r, ht->seqs);
		if (m->status == FRESH && r->contig_locus <= INS_SIZE + SD_INS_SIZE * GRACE_TIMES) {
			g_ptr_array_add(reads, m);
			//p_query("FRESH", m);
		}
	}
	counters = fresh_reads_by_kmer(ht->seqs, reads, ht->o->k);
	tpl_hash::iterator it;
	for (i = 0; i < counters->len; i++) {
		c = (occ_counter*) g_ptr_array_index(counters, i);
		r = &ht->seqs[c->kmer];
		//p_query("START", r);
		b = ext_a_read(ht, all_tpls, t, r, c->count);
		if (b && b->alive) {
			right_tpl_id = right_tpl_to_merge(ht->seqs, b, SM_SIMILARY);
			show_debug_msg(__func__, "RIHGT_TPL_ID: %d \n", right_tpl_id);
			if (right_tpl_id == t->id) {
				//p_tpl_reads(b);
				//p_tpl_reads(t);
				connected = merged_jumped(ht, b, t, N_MISMATCHES);
				break;
			}
		} else {
			rm_global_tpl(all_tpls, t, FRESH);
		}
		free(c);
	}
	g_ptr_array_free(counters, TRUE);
	return connected;
}

void ext_remaining_reads(kmer_t_meta *params) {
	hash_table *ht = params->ht;
	tpl_hash *all_tpls = params->all_tpls;
	int i = 0;
	tpl *t = NULL;
	GPtrArray *tpls_by_pair_pc = sort_tpls_by_pair_pc(params);
	for (i = 0; i < tpls_by_pair_pc->len; i++) {
		t = (tpl*) g_ptr_array_index(tpls_by_pair_pc, i);
		if (!t->alive || t->pair_pc >= 1.0) continue;
		connect_missing_ends(ht, all_tpls, t);
	}
	g_ptr_array_free(tpls_by_pair_pc, TRUE);
}

void kmer_threads(kmer_t_meta *params) {
	GThreadPool *thread_pool = NULL;
	uint64_t i = 0;
	hash_table *ht = params->ht;
	bwa_seq_t *r = NULL, *seqs = ht->seqs;
	occ_counter *counter = NULL;
	GPtrArray *starting_reads = g_ptr_array_sized_new(ht->n_seqs);
	tpl_hash *all_tpls = params->all_tpls;
	tpl *t = NULL;

	show_msg(__func__, "Getting read frequencies ... \n");
	for (i = 0; i < ht->n_seqs; i++) {
		r = &seqs[i];
		//show_debug_msg(__func__, "Query %s: %d\n", r->name, ht->n_kmers[i]);
		if (r->status == FRESH) {
			if (ht->n_kmers[i] > 1) {
				counter = (occ_counter*) malloc(sizeof(occ_counter));
				counter->kmer = i;
				counter->count = ht->n_kmers[i];
				g_ptr_array_add(starting_reads, counter);
			}
		}
	}

	TEST = &seqs[4663319];
	if (!TESTING) {
		show_msg(__func__, "Sorting %d initial reads ... \n", starting_reads->len);
		g_ptr_array_sort(starting_reads, (GCompareFunc) cmp_kmers_by_count);
	}
	show_msg( __func__, "----------- Stage 1: templates without branching ----------\n");
	params->to_try_connect = 1;
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, (gpointer) params, 1, TRUE, NULL);
	for (i = 0; i < starting_reads->len; i++) {
		if (i % 100000 == 0) {
			show_msg(__func__, "%d templates are obtained. \n", params->all_tpls->size());
			show_msg(__func__, "Extending %" ID64 "-th read ... \n", i);
		}
		counter = (occ_counter*) g_ptr_array_index(starting_reads, i);
		//g_thread_pool_push(thread_pool, (gpointer) counter, NULL);
		kmer_ext_thread(counter, params);
		free(counter);
		if (TESTING && fresh_trial >= MAX_FRESH_TRIAL) break;
		//if (i >= 1000000) break;
		//if (all_tpls->size() >= 2) break;
	}
	g_ptr_array_free(starting_reads, TRUE);
	show_msg(__func__, "%d templates are obtained. \n", params->all_tpls->size());

	show_msg( __func__, "----------- Stage 2: connect existing templates ----------\n");
	connect_tpls(params);
	//ext_remaining_reads(params);

	g_thread_pool_free(thread_pool, 0, 1);
}

void save_read_status(hash_table *ht) {
	bwa_seq_t *r = NULL;
	int i = 0;
	char buf[BUFSIZE];
	char *fn = get_output_file("reads.status", kmer_out);
	FILE *f = xopen(fn, "w");
	for (i = 0; i < ht->n_seqs; i++) {
		r = &ht->seqs[i];
		sprintf(buf, "%s\t%d\t%d\n", r->name, r->status, r->contig_id);
		fputs(buf, f);
	}
	free(fn);
	fclose(f);
}

void ext_by_kmers_core(char *lib_file) {
	FILE *contigs = NULL;
	kmer_t_meta *params = (kmer_t_meta*) calloc(1, sizeof(kmer_t_meta));
	tpl_hash all_tpls;
	GPtrArray *read_tpls = NULL, *branching_events = NULL;
	hash_table *ht = NULL;
	char *fn = NULL, *lib_fn_copy = NULL;
	hang_reads = g_ptr_array_sized_new(1024);
	kmer_hash tpl_kmer_hash;

	show_msg(__func__, "Library: %s \n", lib_file);
	lib_fn_copy = str_dup(lib_file);
	ht = load_k_hash(lib_file);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading read hash table: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	branching_events = g_ptr_array_sized_new(BUFSIZ);
	params->ht = ht;
	params->all_tpls = &all_tpls;

	mark_low_qua_reads(ht->seqs, ht->n_seqs);
	kmer_threads(params);

	show_msg(__func__, "Saving contigs: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
	fn = get_output_file("paired.fa", kmer_out);
	contigs = xopen(fn, "w");
	read_tpls = hash_to_array(&all_tpls);
	g_ptr_array_sort(read_tpls, (GCompareFunc) cmp_tpl_by_rev_pair_pc);
	all_tpls.clear();

	save_tpls(read_tpls, contigs, 0, 0, 0);
	save_read_status(ht);
	fflush(contigs);
	fclose(contigs);
	free(fn);

	g_ptr_array_sort(read_tpls, (GCompareFunc) cmp_tpl_by_id);
	get_junction_arr(read_tpls, branching_events);
	g_ptr_array_sort(branching_events, (GCompareFunc) cmp_junc_by_id);
	fn = get_output_file("paired.junctions", kmer_out);
	store_features(fn, branching_events, read_tpls);
	free(fn);

	show_msg(__func__, "Reloading the hash table ... \n");
	reload_table(ht, lib_fn_copy);
	free(lib_fn_copy);
	process_graph(read_tpls, branching_events, ht, kmer_out);
	destroy_ht(ht);
}

int pe_kmer(int argc, char *argv[]) {
	//	validate_junctions("../SRR097897_branch/paired.junctions",
	//			"../SRR097897_branch/paired.fa", "../SRR097897_branch/paired.fa.psl",
	//			"/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897.part.fa");
	//	blat_ref("../SRR097897_part/paired.joint.fa", "../SRR097897_part/paired.joint.fa.psl");
	//	return 0;
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
	if (optind > argc) {
		show_msg(__func__, "Parameters error! \n");
		return 1;
	}
	if (!g_thread_supported())
		g_thread_init( NULL);
	kmer_id_mutex = g_mutex_new();
	tpl_kmers_hash_size <<= kmer_len * 2;

	ext_by_kmers_core(argv[optind]);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));
	return 0;
}
