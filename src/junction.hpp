/*
 * junction.hpp
 *
 *  Created on: 09-May-2013
 *      Author: carl
 */

#ifndef JUNCTION_HPP_
#define JUNCTION_HPP_

#include <glib.h>
#include "hash.cpp"
#include "bwtaln.h"
#include "edge.h"

junction *new_junction(edge *main_tpl, edge *branch_tpl, uint64_t kmer, int locus,
		int ori, int weight);
int find_junc_reads(hash_map *hm, bwa_seq_t *left, bwa_seq_t *right,
		const int max_len, int *weight);
int find_junc_reads_w_tails(hash_map *hm, edge *left, edge *right,
		const int r_shift, const int max_len, int *weight);
void upd_tpl_jun_locus(edge *eg, GPtrArray *branching_events,
		const int kmer_len);
void store_junctions(char *name, GPtrArray *branching_events);

#endif /* JUNCTION_HPP_ */
