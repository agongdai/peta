/*
 * test.hpp
 *
 *  Created on: Jan 3, 2014
 *      Author: carl
 */

#ifndef TEST_HPP_
#define TEST_HPP_

#include <unordered_map>
#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "k_hash.h"
#include "pool.hpp"

using namespace std;

typedef struct {
	int n_match;
	int n_mismatch;
	int n_rep;
	int n_n;
	int n_query_gap;
	int n_query_gap_base;
	int n_ref_gap;
	int n_ref_gap_base;
	char strand;
	char *query;
	int q_len;
	int q_start;
	int q_end;
	char *ref;
	int r_len;
	int r_start;
	int r_end;
	int n_block;
	int *block_size;
	int *query_block_start;
	int *ref_block_start;
} blat_hit;

void test_print_msg(const char *header, const char *fmt, ...);
void test_print_tag(char *msg);
int test_check_missing(tpl *t, int suffix);
int test_to_check(tpl *t, int suffix);
void test_steps(tpl *t, int suffix, char *step);
void test_upd_cursor(tpl *t, int ori);
void test_pool_bad_read_dominates(bwa_seq_t *seqs, tpl *t, pool *p,
		GPtrArray *near_tpls, int ref_c, int c);
int test_check_next_char(tpl *t, int c, int ori);
void test_blat_starting_read(tpl *t, bwa_seq_t *r, int ori);
blat_hit *test_blat_ctg(tpl *t);
void test_cmp_hits(blat_hit *pre_h, blat_hit *next_h);
int test_get_next_ref_char(tpl *t, int ori);
void test_init();
GPtrArray *read_ids(char *ids_fn);
GPtrArray *read_blat_hits(char *blat_psl);
void test_print_blat_hit(blat_hit *h);
void blat_ref(char *joint_fa, char *joint_psl);
void test_tpl_pairs(bwa_seq_t *seqs, tpl *t);
void validate_junctions(char *junc_fn, char *pair_fa, char *pair_psl,
		char *hash_fn);
void read_juncs_from_file(char *junc_fn, char *pair_fa, GPtrArray *all_tpls,
		GPtrArray *all_junctions);
int test_align_tpl_seq(tpl *t, int suffix);

#endif /* TEST_HPP_ */
