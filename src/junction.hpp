/*
 * junction.hpp
 *
 *  Created on: 09-May-2013
 *      Author: carl
 */

#ifndef JUNCTION_HPP_
#define JUNCTION_HPP_

#include <glib.h>
#include <stdint.h>
#include "hash.hpp"
#include "bwtaln.h"
#include "tpl.hpp"

typedef unordered_map<int, int> junc_count;

using namespace std;
typedef struct {
	tpl *main_tpl;
	tpl *branch_tpl;
	int locus;
	int weight;
	uint8_t ori;
	bwa_seq_t *connector; // Kmer when branching
	GPtrArray *reads; // Reads at the junction
	int8_t status; // 0 means good status
} junction;

#ifdef __cplusplus
extern "C" {
#endif

	junction *new_junction(tpl *main_tpl, tpl *branch_tpl, bwa_seq_t *connector,
			int locus, int ori, int weight);
	int tpls_have_junction(tpl *left, tpl *right);
	void p_tpl_juncs(tpl *t, GPtrArray *t_juncs);
	int same_main_juncs(junction *jun_1, junction *jun_2);
	GPtrArray *reset_is_root(tpl *t);
	GPtrArray *find_junc_reads(hash_table *ht, bwa_seq_t *left, bwa_seq_t *right,
			const int max_len, int *weight);
	GPtrArray *find_junc_reads_w_tails(hash_table *ht, tpl *left, tpl *right,
			const int shift, const int max_len, const int ori, int *weight);
	void upd_tpl_jun_locus(tpl *t, GPtrArray *branching_events, const int kmer_len);
	void store_features(char *name, GPtrArray *branching_events,
			GPtrArray *all_tpls);
	gint cmp_junc_by_id(gpointer a, gpointer b);
	void rm_junc_w_dead_tpls(GPtrArray *junctions, tpl *t);
	void disable_tpl_junctions(tpl *t);
	void clean_junctions(GPtrArray *read_tpls, GPtrArray *junctions);
	int branch_on_main(tpl *main, tpl *branch, const int pos, const int mismatches,
			const int exist_ori);
	void destroy_junction(junction *j);
	int count_jun_reads(hash_table *ht, junction *jun);
	bwa_seq_t *get_junc_seq(tpl *left, int l_pos, int *left_len, tpl *right,
			int r_pos, int *right_len, int max_len);
	GPtrArray *find_branch_junctions(GPtrArray *all, tpl *branch);
	GPtrArray *get_nearby_tpls(tpl *t, GPtrArray *reads);
	void p_junction(junction *jun);
	void p_junctions(GPtrArray *juns);
	gint cmp_junc_by_branch_id(gpointer a, gpointer b);
	gint cmp_junc_by_locus(gpointer a, gpointer b);
	GPtrArray *tpl_junctions(tpl *t, GPtrArray *all_juncs, int start_index,
			int to_get_main);
	void filter_junctions(GPtrArray *junctions, GPtrArray *tpls, hash_table *ht);
	int vld_junc_by_mates(tpl *main_tpl, tpl *branch_tpl, GPtrArray *junc_reads,
			hash_table *ht, const int con_pos, const int ins_size, const int ori);

#ifdef __cplusplus
}
#endif

#endif /* JUNCTION_HPP_ */
