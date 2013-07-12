/*
 * junction.c
 *
 *  Created on: 09-May-2013
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
#include "tpl.h"
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "junction.hpp"
#include "graph.hpp"

junction *new_junction(tpl *main_tpl, tpl *branch_tpl, uint64_t kmer,
		int locus, int ori, int weight) {
	junction *j = (junction*) malloc(sizeof(junction));
	j->main_tpl = main_tpl;
	j->branch_tpl = branch_tpl;
	j->locus = main_tpl->len;
	j->ori = ori;
	j->locus = locus;
	j->weight = weight;
	j->kmer = kmer;
	j->reads = g_ptr_array_sized_new(0);
	j->status = 0;
	return j;
}

void destroy_junction(junction *j) {
	if (j) {
		if (j->reads)
			g_ptr_array_free(j->reads, TRUE);
		free(j);
	}
}

void p_tpl_juncs(tpl *t, GPtrArray *t_juncs) {
	show_debug_msg(__func__, "==== Junctions of Template %d ====\n", t->id);
	uint32_t i = 0;
	junction *j = NULL;
	vertex *v = NULL;
	for (i = 0; i < t_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(t_juncs, i);
		p_junction(j);
	}
	show_debug_msg(__func__, "==== Vertexes of Template %d ====\n", t->id);
	for (i = 0; i < t->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(t->vertexes, i);
		printf("\tVertex %d: %d \n", v->id, v->len);
	}
}

gint cmp_junc_by_locus(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->locus) - c_b->locus);
}

gint cmp_junc_by_branch_id(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->branch_tpl->id) - c_b->branch_tpl->id);
}

/**
 * Get junctions whose branch tpl is the given tpl branch
 */
GPtrArray *find_branch_junctions(GPtrArray *all, tpl *branch) {
	uint32_t i = 0;
	junction *j = NULL;
	GPtrArray *hits = g_ptr_array_sized_new(2);
	for (i = 0; i < all->len; i++) {
		j = (junction*) g_ptr_array_index(all, i);
		if (j->branch_tpl == branch)
			g_ptr_array_add(hits, j);
	}
	return hits;
}

int find_junc_reads(hash_map *hm, bwa_seq_t *left, bwa_seq_t *right,
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
	//p_query("Left  seq", left);
	//p_query("Right seq", right);
	//p_query("Junction seq", junction_seq);
	reads = kmer_find_reads(junction_seq, hm, N_MISMATCHES, 0);
	n_reads = reads->len;
	//show_debug_msg(__func__, "# of junction reads: %d \n", n_reads);
	*weight = n_reads;
	bwa_free_read_seq(1, junction_seq);
	g_ptr_array_free(reads, TRUE);
	if (n_reads > 0)
		return 1;
	return 0;
}

bwa_seq_t *get_junc_seq(tpl *left, int l_pos, int *left_len, tpl *right,
		int r_pos, int *right_len, int max_len) {
	bwa_seq_t *junc_seq = blank_seq(max_len);
	bwa_seq_t *left_seq = NULL, *right_seq = NULL;
	left_seq = cut_tpl_tail(left, l_pos, max_len / 2, 0);
	right_seq = cut_tpl_tail(right, r_pos, max_len / 2, 1);
	*left_len = left_seq->len;
	*right_len = right_seq->len;
	memcpy(junc_seq->seq, left_seq->seq, left_seq->len);
	memcpy(junc_seq->seq + left_seq->len, right_seq->seq, right_seq->len);
	junc_seq->len = left_seq->len + right_seq->len;
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return junc_seq;
}

int find_junc_reads_w_tails(hash_map *hm, tpl *main_tpl, tpl *branch_tpl,
		const int pos, const int max_len, const int ori, int *weight) {
	bwa_seq_t *left_seq = NULL, *right_seq = NULL;
	tpl *left_eg = branch_tpl, *right_eg = main_tpl;
	int is_valid = 0, l_pos = branch_tpl->len, r_pos = pos;
	if (ori) {
		left_eg = main_tpl;
		l_pos = pos;
		right_eg = branch_tpl;
		r_pos = 0;
	}
	//show_debug_msg(__func__, "left pos: %d; right pos: %d\n", l_pos, r_pos);
	left_seq = cut_tpl_tail(left_eg, l_pos, max_len / 2, 0);
	//p_query("Left  seq", left_eg->ctg);
	//p_query("Left tail", left_seq);
	right_seq = cut_tpl_tail(right_eg, r_pos, max_len / 2, 1);
	//p_query("Right  seq", right_eg->ctg);
	//p_query("Right tail", right_seq);
	is_valid = find_junc_reads(hm, left_seq, right_seq, max_len, weight);
	bwa_free_read_seq(1, left_seq);
	bwa_free_read_seq(1, right_seq);
	return is_valid;
}

/**
 * Update the junction locus for those tpls connected to itself.
 * Because the locus is not correct when the junction is recorded.
 * Example:
 * Firstly extending to the right and connect to itself:
 *           <<<<<
 *          /     \
 *         /       \
 *         |        ^
 * ------------------
 * Locus:  ^ (value 8)
 * But, later it extends to the left:
 *               <<<<<
 *              /     \
 *             /       \
 *             |        ^
 * ====------------------
 * Locus should be 8 + 4
 */
void upd_tpl_jun_locus(tpl *t, GPtrArray *branching_events, const int kmer_len) {
	int i = 0;
	uint64_t query_int = 0, k = 0;
	junction *jun = NULL;
	if (t->len < kmer_len)
		return;
	for (k = 0; k < branching_events->len; k++) {
		jun = (junction*) g_ptr_array_index(branching_events, k);
		if (jun->main_tpl == t && jun->branch_tpl == t) {
			for (i = 0; i < t->len - kmer_len; i++) {
				query_int = get_kmer_int(t->ctg->seq, i, 1, kmer_len);
				if (query_int == jun->kmer) {
					jun->locus = i;
					break;
				}
			}
		}
	}
}

gint cmp_junctions_by_id(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->main_tpl->id) - c_b->main_tpl->id);
}

void store_features(char *name, GPtrArray *branching_events,
		GPtrArray *all_tpls) {
	junction *jun = NULL;
	uint64_t i = 0;
	tpl *t = NULL;
	char entry[BUFSIZE];
	FILE *f = xopen(name, "w");
	sprintf(entry, "Main\tBranch\tLocus\tWeight\tDirection\n");
	fputs(entry, f);
	//	for (i = 0; i < all_tpls->len; i++) {
	//		t = (tpl*) g_ptr_array_index(all_tpls, i);
	//		sprintf(entry, "[%d, %d]\t[%d, %d]\t0\t%d\t-1\n", t->id, t->len, t->id,
	//				t->len, t->kmer_freq);
	//		fputs(entry, f);
	//	}
	for (i = 0; i < branching_events->len; i++) {
		jun = (junction*) g_ptr_array_index(branching_events, i);
		sprintf(entry, "[%d, %d]\t[%d, %d]\t%d\t%d\t%d\n", jun->main_tpl->id,
				jun->main_tpl->len, jun->branch_tpl->id, jun->branch_tpl->len,
				jun->locus, jun->weight, jun->ori);
		fputs(entry, f);
	}
	fclose(f);
}

void p_junction(junction *jun) {
	show_debug_msg(__func__, "[%d, %d]\t[%d, %d]\t%d\t%d\t%d\n",
			jun->main_tpl->id, jun->main_tpl->len, jun->branch_tpl->id,
			jun->branch_tpl->len, jun->locus, jun->weight, jun->ori);
}

/**
 * Remove those junctions with status not equal to 0 from the list
 */
void remove_dead_junctions(GPtrArray *junctions) {
	junction *j = NULL;
	uint32_t i = 0;
	for (i = 0; i < junctions->len; i++) {
		j = (junction*) g_ptr_array_index(junctions, i);
		if (j->status != 0) {
			destroy_junction(j);
			g_ptr_array_remove_index_fast(junctions, i--);
		}
	}
}

/**
 * Remove templates with no connection to other templates or the length is short
 */
void remove_short_isolate_tpls(GPtrArray *tpls, GPtrArray *junctions, const int max_len) {
	uint32_t i = 0, j = 0;

}

/**
 * Remove duplicate junctions and junctions with dead templates
 */
void clean_junctions(GPtrArray *junctions) {
	junction *junc = NULL, *pre = NULL;
	uint32_t i = 0;
	int not_alive = 0;
	tpl_hash dead_tpls;
	tpl *t = NULL;
	for (i = 0; i < junctions->len; i++) {
		not_alive = 0;
		junc = (junction*) g_ptr_array_index(junctions, i);
		if (!junc->main_tpl->alive) {
			dead_tpls[junc->main_tpl->id] = junc->main_tpl;
			not_alive = 1;
		}
		if (!junc->branch_tpl->alive) {
			dead_tpls[junc->branch_tpl->id] = junc->branch_tpl;
			not_alive = 1;
		}
		if (not_alive) {
			junc->status = 1;
			continue;
		}
		if (pre) {
			if (pre->branch_tpl == junc->branch_tpl && pre->main_tpl
					== junc->main_tpl && pre->kmer == junc->kmer && pre->locus
					== junc->locus && pre->ori == junc->ori) {
				junc->status = 1;
				continue;
			}
		}
		pre = junc;
	}
	//	for (tpl_hash::iterator m = dead_tpls.begin(); m != dead_tpls.end(); ++m) {
	//		t = (tpl*) m->second;
	//		destroy_eg(t);
	//	}
	dead_tpls.clear();
}

/**
 * Main:  ================================
 * Shift:                 ^
 * Ori: 0 (to the right)
 * Branch:                 ---------
 * Check:                  =========
 *                         |||||||||
 *                         ---------
 */
int branch_on_main(const bwa_seq_t *main, const bwa_seq_t *branch,
		const int pos, const int mismatches, const int ori) {
	bwa_seq_t *sub = NULL;
	int similar = 0;
	if (ori) {
		if (pos < branch->len)
			return 0;
		sub = new_seq(main, branch->len, pos - branch->len);
	} else {
		if ((main->len - pos) < branch->len)
			return 0;
		sub = new_seq(main, branch->len, pos);
	}
	p_ctg_seq(__func__, sub);
	p_ctg_seq(__func__, branch);
	//similar = similar_seqs(sub, branch, mismatches, 1, MATCH_SCORE,
	//		MISMATCH_SCORE, INDEL_SCORE);
	similar = seq_ol(sub, branch, branch->len, mismatches);
	free_read_seq(sub);
	show_debug_msg(__func__, "Mismatches: %d; similar: %d\n", mismatches,
			similar);
	return similar;
}

/**
 * Get junctions with the template t as main/branch
 * all_juncs are supposed to be ordered in main/branch template id increasingly
 */
GPtrArray *tpl_junctions(tpl *t, GPtrArray *all_juncs, int start_index,
		int to_get_main) {
	uint32_t i = 0;
	junction *j = NULL;
	GPtrArray *juncs = g_ptr_array_sized_new(32);
	for (i = start_index; i < all_juncs->len; i++) {
		j = (junction*) g_ptr_array_index(all_juncs, i);
		if (to_get_main) {
			if (j->main_tpl == t)
				g_ptr_array_add(juncs, j);
			else
				break;
		} else {
			if (j->branch_tpl == t)
				g_ptr_array_add(juncs, j);
			else
				break;
		}
	}
	return juncs;
}

void filter_branches(GPtrArray *junctions, const int read_len) {
	uint32_t i = 0;
	junction *cur = NULL, *pre = NULL;
	bwa_seq_t *main_seq = NULL, *branch_seq = NULL;
	tpl *main_tpl = NULL, *branch_tpl = NULL;
	for (i = 0; i < junctions->len; i++) {
		cur = (junction*) g_ptr_array_index(junctions, i);
		if (cur->status != 0)
			continue;
		if (!pre) {
			pre = (junction*) g_ptr_array_index(junctions, 0);
			continue;
		}
		if (cur->branch_tpl != pre->branch_tpl) {
			pre = cur;
			continue;
		}
		if (abs(cur->locus - pre->locus) <= 2 && cur->branch_tpl->len <= 2) {
			pre->status = 1;
			cur->status = 1;
			p_junction(pre);
			p_junction(cur);
		} else {
			if (abs(cur->locus - pre->locus) < read_len && pre->ori != cur->ori) {
				branch_seq = branch_tpl->ctg;
			}
		}
		pre = cur;
	}
}

/**
 * Filter out some junctions if:
 * 	- template length 0, connect to the same template at same locus
 *  - hanging branch template can be merged to the main template
 */
void filter_junctions(GPtrArray *junctions, hash_map *hm) {
	uint32_t i = 0, j = 0, start_index = 0;
	junction *junc = NULL;
	tpl *t = NULL;
	GPtrArray *main_junctions = NULL;
	clean_junctions(junctions);
	while(start_index < junctions->len - 1) {
		show_debug_msg(__func__, "======== Started at %d =======\n", start_index);
		junc = (junction*) g_ptr_array_index(junctions, start_index);
		t = junc->main_tpl;
		main_junctions = tpl_junctions(t, junctions, start_index, 1);
		filter_branches(main_junctions, hm->o->read_len);
		start_index += main_junctions->len;
		g_ptr_array_free(main_junctions, TRUE);
	}
}
