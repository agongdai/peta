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
#include "edge.h"

using namespace std;
typedef struct {
	edge *main_tpl;
	edge *branch_tpl;
	int locus;
	int weight;
	uint8_t ori;
	uint64_t kmer;	// Kmer when branching
} junction;

#ifdef __cplusplus
extern "C" {
#endif

	junction *new_junction(edge *main_tpl, edge *branch_tpl, uint64_t kmer, int locus,
			int ori, int weight);
	int find_junc_reads(hash_map *hm, bwa_seq_t *left, bwa_seq_t *right,
			const int max_len, int *weight);
	int find_junc_reads_w_tails(hash_map *hm, edge *left, edge *right,
			const int shift, const int max_len, const int ori, int *weight);
	void upd_tpl_jun_locus(edge *eg, GPtrArray *branching_events,
			const int kmer_len);
	void store_junctions(char *name, GPtrArray *branching_events);
	gint cmp_junctions_by_id(gpointer a, gpointer b);

#ifdef __cplusplus
}
#endif

#endif /* JUNCTION_HPP_ */
