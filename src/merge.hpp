/*
 * merge.hpp
 *
 *  Created on: Jul 29, 2013
 *      Author: carl
 */

#ifndef MERGE_HPP_
#define MERGE_HPP_

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <glib.h>
#include "utils.h"
#include "tpl.hpp"
#include "peseq.h"

void mv_reads_bt_tpls(tpl *from, tpl *to, int ol, int ori);
int merge_tpls(tpl *left, tpl *right, int ol, int rev_com);
int merged_jumped(hash_table *ht, tpl *t, tpl *jumped, int mis);
void merge_tpl_to_right(tpl *t, tpl *jumped, int ol, int rev_com);
void merge_tpl_to_left(tpl *t, tpl *jumped, int ol, int rev_com);
int connect_at_locus_right(hash_table *ht, tpl *t, tpl *b, int t_locus, int b_locus);
int connect_both_ends(hash_table *ht, GPtrArray *anchors, tpl *t);
int connect_one_end(hash_table *ht, GPtrArray *anchors, tpl *t);
int left_tpl_is_paired(bwa_seq_t *seqs, tpl *left, tpl *right, float pair_pc);
int right_tpl_to_merge(bwa_seq_t *seqs, tpl *left, float pair_pc);

#endif /* MERGE_HPP_ */
