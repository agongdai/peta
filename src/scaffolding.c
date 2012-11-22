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
 * Return 0: the two edges are in the same direction
 */
int get_edges_ori(edge *eg_left, edge *eg_right) {
	GPtrArray *paired_reads = NULL;
	int n_forward = 0, n_backward = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	int i = 0, ori = 0;
	paired_reads = find_unconditional_paired_reads(eg_left->reads,
			eg_right->reads);
	p_readarray(paired_reads, 1);
	if (paired_reads->len > 0) {
		for (i = 0; i < paired_reads->len - 1; i += 2) {
			show_debug_msg(__func__, "i: %d/%d \n", i, paired_reads->len);
			read = g_ptr_array_index(paired_reads, i);
			mate = g_ptr_array_index(paired_reads, i + 1);
			//		p_query(__func__, read);
			//		p_query(__func__, mate);
			if (read->rev_com == mate->rev_com)
				n_forward++;
			else
				n_backward++;
		}
	}
	if (n_forward > n_backward) {
		eg_right->ori = eg_left->ori;
	} else {
		eg_right->ori = eg_left->ori ? 0 : 1;
		ori = 1;
	}
	g_ptr_array_free(paired_reads, TRUE);
	return ori;
}

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
	int i = 0, forward_count = 0, backward_count = 0, min = 0, ori = 0;
	bwa_seq_t *read = NULL, *mate = NULL;

	ori = get_edges_ori(eg_1, eg_2);
	paired_reads = find_unconditional_paired_reads(eg_1->reads, eg_2->reads);
	//	show_debug_msg(__func__, "Ordering edge %d and %d: %d pairs ...\n", eg_1->id,
	//			eg_2->id, paired_reads->len);
	if (paired_reads->len < MIN_VALID_PAIRS) {
		g_ptr_array_free(paired_reads, TRUE);
		return 0;
	} else {
		orders = (int*) calloc(paired_reads->len / 2 + 1, sizeof(int));
		for (i = 0; i < paired_reads->len; i += 2) {
			read = g_ptr_array_index(paired_reads, i);
			mate = g_ptr_array_index(paired_reads, i + 1);
			if (ori && read->rev_com == mate->rev_com)
				continue;
			if (!ori && read->rev_com != mate->rev_com)
				continue;
			if (is_left_mate(read->name)) {
				orders[i / 2] = (read->rev_com) ? -1 : 1;
			} else {
				orders[i / 2] = (read->rev_com) ? 1 : -1;
			}
		}
	}
	for (i = 0; i < paired_reads->len / 2; i++) {
		if (orders[i] == 1)
			forward_count++;
		else
			backward_count++;
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

GPtrArray *get_probable_in_out(GPtrArray *all_edges, edge *eg, bwa_seq_t *seqs) {
	edgearray *probable_in_out = NULL;
	int i = 0, n_in_out = 0;
	edge *in_out = NULL;
	bwa_seq_t *read = NULL, *mate = NULL, *rev = NULL;
	int *edge_index_got = NULL;
	edge_index_got = (int*) calloc(all_edges->len + 1, sizeof(int));
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		if (read->status == USED)
			continue;
		mate = get_mate(read, seqs);
		if (mate->status == FRESH || mate->contig_id < 0)
			continue;
		//p_query(__func__, mate);
		if (edge_index_got[mate->contig_id] == 0)
			n_in_out++;
		//show_debug_msg(__func__, "n_in_out: %d \n", n_in_out);
		edge_index_got[mate->contig_id] = 1;
	}
	probable_in_out = g_ptr_array_sized_new(n_in_out + 1);
	for (i = 0; i < all_edges->len; i++) {
		if (edge_index_got[i] == 1) {
			in_out = g_ptr_array_index(all_edges, i);
			rev = new_mem_rev_seq(eg->contig, eg->len, 0);
			if (!has_reads_in_common(eg, in_out) && !share_subseq(eg->contig,
					in_out->contig, MISMATCHES, 100) && !share_subseq(rev,
							in_out->contig, MISMATCHES, 100))
				g_ptr_array_add(probable_in_out, in_out);
			bwa_free_read_seq(1, rev);
		}
	}
	free(edge_index_got);
	return probable_in_out;
}

/*
 * Check whether two edges contain any same read
 */
int has_reads_in_common(edge *eg_1, edge *eg_2) {
	int i = 0, j = 0;
	bwa_seq_t *read = NULL, *same = NULL;
	g_ptr_array_sort(eg_1->reads, (GCompareFunc) cmp_read_by_name);
	g_ptr_array_sort(eg_2->reads, (GCompareFunc) cmp_read_by_name);
	for (i = 0; i < eg_1->reads->len; i++) {
		read = g_ptr_array_index(eg_1->reads, i);
		for (j = i; j < eg_2->reads->len; j++) {
			same = g_ptr_array_index(eg_2->reads, j);
			if (read == same)
				return 1;
			if (atoi(read->name) < atoi(same->name))
				break;
		}
	}
	return 0;
}

edge *merge_edges(edge *eg_1, edge *eg_2) {
	int order = 0, ori = 0, i = 0;
	bwa_seq_t *read = NULL;
	order = order_two_edges(eg_1, eg_2);
	ori = get_edges_ori(eg_1, eg_2);
	if (order == 0)
		return 0;
	if (order == 1) {
		if (ori == 1) { // Means two edges are not in the same direction
			seq_reverse(eg_2->len, eg_2->contig->seq, 1);
			for (i = 0; i < eg_2->reads->len; i++) {
				read = g_ptr_array_index(eg_2->reads, i);
				read->rev_com = read->rev_com ? 0 : 1;
			}
		}
		merge_seq_to_left(eg_1->contig, eg_2->contig, 0);
		eg_1->len = eg_1->contig->len;
		concat_readarray(eg_1->reads, eg_2->reads);
		concat_readarray(eg_1->pairs, eg_2->pairs);
		destroy_eg(eg_2);
		upd_reads(eg_1, MISMATCHES);
		return eg_1;
	} else {
		if (ori == 1) { // Means two edges are not in the same direction
			seq_reverse(eg_1->len, eg_1->contig->seq, 1);
			for (i = 0; i < eg_1->reads->len; i++) {
				read = g_ptr_array_index(eg_1->reads, i);
				read->rev_com = read->rev_com ? 0 : 1;
			}
		}
		merge_seq_to_left(eg_2->contig, eg_1->contig, 0);
		eg_2->len = eg_2->contig->len;
		concat_readarray(eg_2->reads, eg_1->reads);
		concat_readarray(eg_2->pairs, eg_1->pairs);
		destroy_eg(eg_1);
		upd_reads(eg_2, MISMATCHES);
		return eg_2;
	}
}

GPtrArray *scaffolding(GPtrArray *single_edges, const int insert_size,
		bwa_seq_t *seqs) {
	edge *eg_i = NULL, *eg_j = NULL;
	int i = 0, j = 0, gap = 0, order = 0;
	edgearray *probable_in_out = NULL;

	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		//show_debug_msg(__func__, "Trying edge [%d/%d, %d] \n", eg_i->id, single_edges->len, eg_i->len);
		probable_in_out = get_probable_in_out(single_edges, eg_i, seqs);
		for (j = 0; j < probable_in_out->len; j++) {
			eg_j = g_ptr_array_index(probable_in_out, j);
			if (eg_i == eg_j)
				continue;
			order = order_two_edges(eg_i, eg_j);
			if (order != 0) {
				gap = est_pair_gap(eg_i, eg_j, order, insert_size);
				gap = 0; // @TODO: orientation of edge is not determined yet, skip gap size estimation
				if (gap >= 0) {
					show_debug_msg(__func__,
							"Order of edge %d and edge %d: %d; Gap: %d\n",
							eg_i->id, eg_j->id, order, gap);
					if (order == 1) {
						g_ptr_array_add(eg_i->out_egs, eg_j);
						g_ptr_array_add(eg_j->in_egs, eg_i);
					} else {
						g_ptr_array_add(eg_j->in_egs, eg_i);
						g_ptr_array_add(eg_i->out_egs, eg_j);
					}
				}
			}
		}
		g_ptr_array_free(probable_in_out, TRUE);
	}
	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		if (eg_i && eg_i->in_egs->len == 0 && eg_i->alive) {
			eg_i->is_root = 1;
		}
		p_flat_eg(eg_i);
	}
	return NULL;
}
