/*
 * graph.cpp
 *
 *  Created on: Jun 30, 2013
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
#include "kmers.hpp"
#include "ass.hpp"
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "junction.hpp"
#include "hash.hpp"

void assign_reads2tpl(edge *eg, hash_map *hm) {
	int j = 0, i = 0, read_len = hm->o->read_len;
	bwa_seq_t *seq = NULL, *read = NULL;
	GPtrArray *hits = NULL;

	for (i = 0; i <= eg->len - read_len; i++) {
		seq = new_seq(eg->ctg, read_len, i);
		hits = align_full_seq(seq, hm, 2);
		for (j = 0; j < hits->len; j++) {
			read = (bwa_seq_t*) g_ptr_array_index(hits, j);
			add_read2edge(eg, read);
		}
		bwa_free_read_seq(1, seq);
		g_ptr_array_free(hits, TRUE);
	}
}

gint cmp_junc_by_locus(gpointer a, gpointer b) {
	junction *c_a = *((junction**) a);
	junction *c_b = *((junction**) b);
	return ((c_a->locus) - c_b->locus);
}

void break_tpls(GPtrArray *all_tpls, GPtrArray *junctions) {
	edge *eg = NULL;
	junction *junc = NULL;
	uint64_t i = 0, j = 0;
	GPtrArray *eg_junc = NULL;
	for (i = 0; i < all_tpls->len; i++) {
		eg = (edge*) g_ptr_array_index(all_tpls, i);
		eg_junc = g_ptr_array_sized_new(4);
		for (j = 0; j < junctions->len; j++) {
			junc = (junction*) g_ptr_array_index(junctions, j);
			if (junc->main_tpl == eg)
				g_ptr_array_add(eg_junc, junc);
		}
		g_ptr_array_sort(eg_junc, (GCompareFunc) cmp_junc_by_locus);
		for (j = 0; j < eg_junc->len; j++) {
			junc = (junction*) g_ptr_array_index(eg_junc, j);

		}
		g_ptr_array_free(eg_junc, TRUE);
	}
}
