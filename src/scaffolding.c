/*
 * scaffolding.c
 *
 *  Created on: 18-Nov-2012
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "rand.h"
#include "ass.h"
#include "bwase.h"
#include "utils.h"
#include "peseq.h"
#include "pool.h"
#include "pechar.h"
#include "pealn.h"
#include "roadmap.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"
#include "pelib.h"
#include "pepath.h"
#include "scaffolding.h"

/**
 * order:
 * 	1: eg_1 is left
 *  -1: eg_1 is right
 */
int est_pair_gap(edge *eg_1, edge *eg_2, const int order, const int insert_size) {
	GPtrArray *paired_reads = NULL;
	int i = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	int len = 0, gap_sum = 0, gap_ave = -1, n_counted = 0;
	paired_reads = find_unconditional_paired_reads(eg_1->reads, eg_2->reads);
	for (i = 0; i < paired_reads->len; i += 2) {
		read = g_ptr_array_index(paired_reads, i);
		mate = g_ptr_array_index(paired_reads, i + 1);
		if (read->rev_com != mate->rev_com)
			continue;
		n_counted++;
		if (order == 1) {
			len = (eg_1->len - read->shift) + mate->shift;
		} else
			len = (eg_2->len - mate->shift) + read->shift;
		gap_sum += insert_size - len;
	}
	if (n_counted > 0) {
		gap_ave = gap_sum / n_counted;
		if (gap_ave < 0 - insert_size / 2 || gap_ave > insert_size) {
			gap_ave = -1;
		} else {
			for (i = 0; i < paired_reads->len; i++) {
				read = g_ptr_array_index(paired_reads, i);
				read->status = USED; // Mark this read as used, not used to connect edges anymore.
			}
			gap_ave = gap_ave > 0 ? gap_ave : 0;
		}
	}
	g_ptr_array_free(paired_reads, TRUE);
	return gap_ave;
}

/**
 * Return:
 * 	0: They cannot be put together
 * 	1: eg_1 is left, eg_2 is right
 *  -1: eg_1 is right, eg_2 is left
 */
int order_two_edges(edge *eg_1, edge *eg_2) {
	GPtrArray *paired_reads = NULL;
	int *orders = NULL, consensus_order = 0;
	int i = 0, forward_count = 0, backward_count = 0, min = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	paired_reads = find_unconditional_paired_reads(eg_1->reads, eg_2->reads);
	show_debug_msg(__func__, "Ordering edge %d and %d: %d pairs ...\n", eg_1->id,
			eg_2->id, paired_reads->len);
	if (paired_reads->len < MIN_VALID_PAIRS)
		orders = 0;
	else {
		orders = (int*) calloc(paired_reads->len / 2 + 1, sizeof(int));
		for (i = 0; i < paired_reads->len; i += 2) {
			read = g_ptr_array_index(paired_reads, i);
			mate = g_ptr_array_index(paired_reads, i + 1);
			if (read->rev_com != mate->rev_com)
				continue;
			if (is_left_mate(read->name)) {
				orders[i / 2] = (read->rev_com) ? -1 : 1;
			} else {
				orders[i / 2] = (read->rev_com) ? 1 : -1;
			}
			break;
		}
		for (i = 0; i < paired_reads->len / 2; i++) {
			if (orders[i] == 1)
				forward_count++;
			else
				backward_count++;
		}
	}
	g_ptr_array_free(paired_reads, TRUE);
	free(orders);
	min = forward_count > backward_count ? backward_count : forward_count;
	consensus_order = forward_count > backward_count ? 1 : -1;
	if (min <= 2) {
		return consensus_order;
	}
	return 0;
}

GPtrArray *scaffolding(GPtrArray *single_edges, const int insert_size) {
	edge *eg_i = NULL, *eg_j = NULL;
	int i = 0, j = 0, gap = 0, order = 0;

	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		for (j = 0; j < single_edges->len; j++) {
			eg_j = g_ptr_array_index(single_edges, j);
			if (eg_i == eg_j)
				continue;
			order = order_two_edges(eg_i, eg_j);
			if (order != 0) {
				gap = est_pair_gap(eg_i, eg_j, order, insert_size);
				gap = 1; // @TODO: temp
				if (gap >= 0) {
					show_debug_msg(__func__,
							"Order of edge %d and edge %d: %d; Gap: %d\n",
							eg_i->id, eg_j->id, order, gap);
					if (order == 1) {
						g_ptr_array_add(eg_i->out_egs, eg_j);
					} else
						g_ptr_array_add(eg_j->in_egs, eg_i);
				}
			}
		}
	}
}
