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
int get_edges_ori(edge *eg_left, edge *eg_right, bwa_seq_t *seqs) {
	GPtrArray *paired_reads = NULL;
	int n_forward = 0, n_backward = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	int i = 0, ori = 0;
	paired_reads = find_unconditional_paired_reads(eg_left,
			eg_right, seqs);
	if (paired_reads->len > 0) {
		for (i = 0; i < paired_reads->len - 1; i += 2) {
			read = g_ptr_array_index(paired_reads, i);
			mate = g_ptr_array_index(paired_reads, i + 1);
			if (read->rev_com == mate->rev_com)
				n_forward++;
			else
				n_backward++;
		}
	}
	if (n_forward > 0 || n_backward > 0) {
		if (n_forward > n_backward) {
			eg_right->ori = eg_left->ori;
		} else {
			eg_right->ori = eg_left->ori ? 0 : 1;
			ori = 1;
		}
	}
	g_ptr_array_free(paired_reads, TRUE);
	return ori;
}

/**
 * order:
 * 	1: eg_1 is left
 *  -1: eg_1 is right
 */
int est_pair_gap(edge *eg_1, edge *eg_2, const int order, const int insert_size, bwa_seq_t *seqs) {
	GPtrArray *paired_reads = NULL;
	int i = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	int len = 0, gap_sum = 0, gap_ave = -1, n_counted = 0;
	paired_reads = find_unconditional_paired_reads(eg_1, eg_2, seqs);
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
int order_two_edges(edge *eg_1, edge *eg_2, bwa_seq_t *seqs) {
	GPtrArray *paired_reads = NULL;
	int *orders = NULL, consensus_order = 0;
	int i = 0, forward_count = 0, backward_count = 0, min = 0, ori = 0;
	bwa_seq_t *read = NULL, *mate = NULL;

	ori = get_edges_ori(eg_1, eg_2, seqs);
	paired_reads = find_unconditional_paired_reads(eg_1, eg_2, seqs);
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
		if (mate->status == FRESH || mate->contig_id < 0 || read->contig_id == mate->contig_id)
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
	int i = 0, exists = 0;
	bwa_seq_t *read = NULL;
	// Iterate through the array with fewer reads, to perform binary search in the array with more reads.
	edge *eg_fewer = eg_1, *eg_more = eg_2;
	if (eg_1->len > eg_2->len) {
		eg_fewer = eg_2;
		eg_more = eg_1;
		g_ptr_array_sort(eg_1->reads, (GCompareFunc) cmp_read_by_name);
	} else {
		g_ptr_array_sort(eg_2->reads, (GCompareFunc) cmp_read_by_name);
	}
	for (i = 0; i < eg_fewer->reads->len; i++) {
		read = g_ptr_array_index(eg_fewer->reads, i);
		exists = binary_exists(eg_more->reads, read);
		if (exists) {
			//show_debug_msg(__func__, "Edge [%d, %d] and [%d, %d] share a common read. \n", eg_1->id, eg_1->len, eg_2->id, eg_2->len);
			p_query(__func__, read);
			return 1;
		}
	}
	return 0;
}

edge *merge_two_ol_edges(hash_table *ht, edge *eg_1, edge *eg_2, const int ol) {
	int i = 0;
	bwa_seq_t *r = NULL;

	//show_debug_msg(__func__, "Merging edges [%d, %d] and [%d, %d] by overlapping... \n", eg_1->id, eg_1->len, eg_2->id, eg_2->len);
	eg_1->len -= ol;
	eg_1->contig->len = eg_1->len;
	eg_1->contig->seq[eg_1->len] = '\0';
	merge_seq_to_left(eg_1->contig, eg_2->contig, 0);
	eg_1->len = eg_1->contig->len;
	//show_debug_msg(__func__, "Concating reads ... \n");
	concat_readarray(eg_1->reads, eg_2->reads);
	//show_debug_msg(__func__, "Concating pairs ... \n");
	concat_readarray(eg_1->pairs, eg_2->pairs);
	//show_debug_msg(__func__, "Clearing reads ... \n");
	clear_used_reads(eg_2, 0);
	//show_debug_msg(__func__, "Updating reads %d=>%d ... \n", eg_1->id, eg_1->reads->len);
	upd_reads_by_ht(ht, eg_1, MISMATCHES);
	//show_debug_msg(__func__, "Updated ... \n");
	for (i = 0; i < eg_1->reads->len; i++) {
		r = g_ptr_array_index(eg_1->reads, i);
		r->status = TRIED;
		r->contig_id = eg_1->id;
	}
	for (i = 0; i < eg_1->pairs->len; i++) {
		r = g_ptr_array_index(eg_1->pairs, i);
		r->status = USED;
	}
	eg_2->alive = 0;
	return eg_1;
}

edge *merge_edges(edge *eg_1, edge *eg_2, hash_table *ht) {
	int order = 0, ori = 0, i = 0;
	bwa_seq_t *read = NULL;
	order = order_two_edges(eg_1, eg_2, ht->seqs);
	ori = get_edges_ori(eg_1, eg_2, ht->seqs);
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
		merge_two_ol_edges(ht, eg_1, eg_2, 0);
		destroy_eg(eg_2);
		return eg_1;
	} else {
		if (ori == 1) { // Means two edges are not in the same direction
			seq_reverse(eg_1->len, eg_1->contig->seq, 1);
			for (i = 0; i < eg_1->reads->len; i++) {
				read = g_ptr_array_index(eg_1->reads, i);
				read->rev_com = read->rev_com ? 0 : 1;
			}
		}
		merge_two_ol_edges(ht, eg_2, eg_1, 0);
		destroy_eg(eg_1);
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
			order = order_two_edges(eg_i, eg_j, seqs);
			if (order != 0) {
				//gap = est_pair_gap(eg_i, eg_j, order, insert_size);
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
	}
	return NULL;
}

/**
 * Merge the edges with some overlapping
 */
void merge_ol_edges(edgearray *single_edges, const int insert_size,
		hash_table *ht) {
	edge *eg_i = NULL, *eg_j = NULL;
	int i = 0, j = 0, some_one_merged = 1, ol = 0, has_common_read = 0;
	bwa_seq_t *seqs = ht->seqs, *rev = NULL;
	readarray *paired_reads = NULL;
	int rl = 0;
	clock_t t = clock();

	rl = seqs->len;

	while (some_one_merged) {
		some_one_merged = 0;
		for (i = 0; i < single_edges->len; i++) {
			eg_i = g_ptr_array_index(single_edges, i);
			eg_i->visited = 1;
			if (!eg_i->alive)
				continue;
			show_debug_msg(__func__, "Trying edge [%d/%d, %d]... %.2f sec\n", eg_i->id, single_edges->len, eg_i->len, (float) (clock() - t) / CLOCKS_PER_SEC);
			for (j = 0; j < single_edges->len; j++) {
				eg_j = g_ptr_array_index(single_edges, j);
				eg_i->visited = 1;
				if (!eg_j->alive || eg_i == eg_j || (eg_i->len < 100
						&& eg_j->len < 100))
					continue;
				if (eg_i->visited && eg_j->visited)
					continue;
				paired_reads = find_unconditional_paired_reads(eg_i,
						eg_j, seqs);
				if (paired_reads->len == 0) {
					g_ptr_array_free(paired_reads, TRUE);
					continue;
				}
				g_ptr_array_free(paired_reads, TRUE);
				rev = new_mem_rev_seq(eg_j->contig, eg_j->contig->len, 0);

				//show_debug_msg(__func__, "\t sub edge [%d/%d, %d]... %.2f sec\n", eg_j->id, single_edges->len, eg_j->len, (float) (clock() - t) / CLOCKS_PER_SEC);
				ol = find_ol(eg_i->contig, eg_j->contig, MISMATCHES);
				// 1. If the overlapping length is shorter than read length,
				// 		We expect that no common reads
				// 2. If the overlapping length is longer than read length,
				//		There mush be some common reads.
				has_common_read = has_reads_in_common(eg_i, eg_j); 
				if ((ol >= MATE_OVERLAP_THRE && ol < rl
						&& !has_common_read) || (ol > rl
						&& has_common_read)) {
					show_debug_msg(__func__,
							"Merging edge [%d, %d] to edge [%d, %d]\n",
							eg_j->id, eg_j->len, eg_i->id, eg_i->len);
					merge_two_ol_edges(ht, eg_i, eg_j, ol);
					eg_i->visited = 0;
					some_one_merged = 1;
				} else {
					ol = find_ol(eg_i->contig, rev, MISMATCHES);
					if (ol >= MATE_OVERLAP_THRE) {
						show_debug_msg(__func__,
								"Merging edge [%d, %d] to edge [%d, %d]\n",
								eg_j->id, eg_j->len, eg_i->id, eg_i->len);
						bwa_free_read_seq(1, eg_j->contig);
						eg_j->contig = rev;
						merge_two_ol_edges(ht, eg_i, eg_j, ol);
						eg_i->visited = 0;
						some_one_merged = 1;
						continue; // In case the 'rev' is freed accidently.
					}
				}
				bwa_free_read_seq(1, rev);
			}
		}
	}
	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		eg_i->visited = 0;
		if (!eg_i->alive) {
			destroy_eg(eg_i);
			g_ptr_array_remove_index_fast(single_edges, i);
			i--;
		}
	}
}
