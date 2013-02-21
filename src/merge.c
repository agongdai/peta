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
#include "merge.h"
#include "pehash.h"
#include "hits.h"

GMutex *edge_mutex = NULL;

/**
 * Return 0: the two edges are in the same direction
 */
int get_edges_ori(edge *eg_left, edge *eg_right, bwa_seq_t *seqs) {
	GPtrArray *paired_reads = NULL;
	int n_forward = 0, n_backward = 0;
	bwa_seq_t *read = NULL, *mate = NULL;
	int i = 0, ori = 0;
	paired_reads = find_unconditional_paired_reads(eg_left, eg_right, seqs);
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
	//show_debug_msg(__func__, "Ordering edge %d and %d: %d pairs ...\n",
	//		eg_1->id, eg_2->id, paired_reads->len);
	// p_readarray(paired_reads, 1);
	if (paired_reads->len < MIN_VALID_PAIRS) {
		g_ptr_array_free(paired_reads, TRUE);
		return 0;
	} else {
		orders = (int*) calloc(paired_reads->len / 2 + 2, sizeof(int));
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
	}
	for (i = 0; i < eg_fewer->reads->len; i++) {
		read = g_ptr_array_index(eg_fewer->reads, i);
		exists = binary_exists(eg_more->reads, read);
		if (exists) {
			return 1;
		}
	}
	return 0;
}

GPtrArray *get_probable_in_out(GPtrArray *all_edges, edge *eg, bwa_seq_t *seqs) {
	edgearray *probable_in_out = NULL, *raw_in_outs = NULL;
	int i = 0;
	edge *in_out = NULL;
	bwa_seq_t *read = NULL, *mate = NULL;
	raw_in_outs = g_ptr_array_sized_new(all_edges->len + 1);
	//show_debug_msg(__func__,
	//		"Checking probable in_out edges of edge [%d, %d]...\n", eg->id,
	//		eg->len);
	for (i = 0; i < eg->reads->len; i++) {
		read = g_ptr_array_index(eg->reads, i);
		if (read->status == USED)
			continue;
		mate = get_mate(read, seqs);
		if (mate->status == FRESH || mate->contig_id < 0 || read->contig_id
				== mate->contig_id)
			continue;
		//show_debug_msg(__func__, "mate contig id: %d/%d \n", mate->contig_id,
		//		all_edges->len);
		in_out = edgearray_find_id(all_edges, mate->contig_id);
		if (in_out) {
			g_ptr_array_uni_add(raw_in_outs, in_out);
		}
	}
	probable_in_out = g_ptr_array_sized_new(raw_in_outs->len + 1);
	//show_debug_msg(__func__, "Checking shared subseq... \n");
	for (i = 0; i < raw_in_outs->len; i++) {
		in_out = g_ptr_array_index(raw_in_outs, i);
		//show_debug_msg(__func__, "[%d] Probable in_out: [%d, %d] \n", i,
		//		in_out->id, in_out->len);
		// p_ctg_seq("Contig", eg->contig);
		if (!share_subseq_byte(eg->contig->seq, eg->len, in_out->contig,
				MISMATCHES, 100) && !share_subseq_byte(eg->contig->rseq,
				eg->len, in_out->contig, MISMATCHES, 100)
				&& !has_reads_in_common(eg, in_out)) {
			in_out->tid = eg->tid;
			g_ptr_array_add(probable_in_out, in_out);
		}
	}
	g_ptr_array_free(raw_in_outs, TRUE);
	return probable_in_out;
}

/**
 * Merge two overlapped edges to a single one
 */
edge *merge_two_ol_edges(reads_ht *rht, hash_table *ht, edge *eg_1, edge *eg_2,
		const int ol) {
	int i = 0;
	bwa_seq_t *r = NULL, *mate = NULL;

	if (!eg_1->alive || !eg_2->alive)
		return eg_1;

	show_debug_msg(__func__,
			"Merging edges [%d, %d] and [%d, %d] by overlapping... \n",
			eg_1->id, eg_1->len, eg_2->id, eg_2->len);
	eg_1->len -= ol;
	eg_1->contig->len = eg_1->len;
	eg_1->contig->seq[eg_1->len] = '\0';
	rm_read_from_ht(rht, eg_1->contig);
	rm_read_from_ht(rht, eg_2->contig);
	merge_seq_to_left(eg_1->contig, eg_2->contig, 0);
	//show_debug_msg(__func__, "Concating reads ... \n");
	for (i = 0; i < eg_2->reads->len; i++) {
		r = g_ptr_array_index(eg_2->reads, i);
		if (r->contig_id != eg_1->id) {
			g_ptr_array_add(eg_1->reads, r);
			r->contig_id = eg_1->id;
			r->shift += eg_2->len - ol;
		}
	}
	eg_1->len = eg_1->contig->len;
	//show_debug_msg(__func__, "Clearing reads ... \n");
	clear_used_reads(eg_2, 0);
	//show_debug_msg(__func__, "Updated ... \n");
	for (i = 0; i < eg_1->reads->len; i++) {
		r = g_ptr_array_index(eg_1->reads, i);
		r->status = TRIED;
	}
	for (i = 0; i < eg_1->reads->len; i++) {
		r = g_ptr_array_index(eg_1->reads, i);
		mate = get_mate(r, ht->seqs);
		if (mate->contig_id == r->contig_id && mate->status != USED
				&& r->status != USED) {
			g_ptr_array_add(eg_1->pairs, r);
			g_ptr_array_add(eg_1->pairs, mate);
			r->status = USED;
			mate->status = USED;
		}
	}
	eg_2->alive = 0;
	add_read_to_ht(rht, eg_1->contig);
	set_rev_com(eg_1->contig);
	g_ptr_array_sort(eg_1->reads, (GCompareFunc) cmp_read_by_name);
	return eg_1;
}

typedef struct {
	edgearray *single_edges;
	hash_table *ht;
	int insert_size;
	int sd_insert_size;
	int start;
	int end;
	int tid;
	reads_ht *rht;
} merging_paras_t;

int check_insert_size(edge *eg_left, edge *eg_right, readarray *paired_reads,
		const int insert_size, const int sd_insert_size, const int ol,
		const int type) {
	int i = 0, dis = 0;
	bwa_seq_t *read_left = NULL, *read_right = NULL, *tmp = NULL;
	for (i = 0; i < paired_reads->len - 1; i += 2) {
		read_left = g_ptr_array_index(paired_reads, i);
		read_right = g_ptr_array_index(paired_reads, i + 1);
		if (read_left->contig_id != eg_left->id) {
			tmp = read_left;
			read_left = read_right;
			read_right = tmp;
		}
		// Tail-Head
		if (type == TAIL_HEAD) {
			dis = (eg_left->len - read_left->shift) + (read_right->shift) - ol;
		}
		// Tail-Tail
		if (type == TAIL_TAIL) {
			dis = (eg_left->len - read_left->shift) + (read_right->len
					- read_right->shift) - ol;
		}
		// Head-Head
		if (type == HEAD_HEAD) {
			dis = (read_left->shift) + (read_right->shift) - ol;
		}
		if (dis > insert_size + SD_TIMES * sd_insert_size)
			return 0;
	}
	return 1;
}

/**
 * For each pair of contigs, check:
 * 1: -------------->
 * 2:			--------------->
 *
 * 1:             -------------->
 * 2: --------------->
 *
 * 1: <--------------
 * 2: 			--------------->
 *
 * In next round, edge 1 and 2 maybe exchanged
 */
void *merge_ol_edges_thread(void *data) {
	edge *eg_i = NULL, *eg_j = NULL;
	int i = 0, j = 0, some_one_merged = 1, ol = 0, nm = 0, to_merge = 0;
	bwa_seq_t *seqs = NULL, *rev = NULL;
	readarray *paired_reads = NULL;
	GPtrArray *edge_candidates = NULL;
	merging_paras_t *d = (merging_paras_t*) data;
	hash_opt *hash_o = d->ht->o;

	seqs = d->ht->seqs;

	while (some_one_merged) {
		some_one_merged = 0;
		for (i = d->start; i < d->end; i++) {
			eg_i = g_ptr_array_index(d->single_edges, i);
			eg_i->tid = d->tid;
			if (!eg_i->alive)
				continue;
			if (edge_candidates) {
				g_ptr_array_free(edge_candidates, TRUE);
				edge_candidates = NULL;
			}
			edge_candidates = find_edges_ol(d->rht, eg_i->contig,
					d->single_edges);
			for (j = 0; j < edge_candidates->len; j++) {
				bwa_free_read_seq(1, rev);
				rev = NULL;
				to_merge = 0;
				eg_j = g_ptr_array_index(edge_candidates, j);
				if (!eg_j->alive || eg_i == eg_j || (eg_i->len < 100
						&& eg_j->len < 100))
					continue;
				if (eg_i->visited && eg_j->visited)
					continue;

				ol = find_ol(eg_i->contig, eg_j->contig, MAX_EDGE_NM);
				nm = get_mismatches_on_ol(eg_i->contig, eg_j->contig, ol,
						MAX_EDGE_NM);
				/**
				 p_ctg_seq(__func__, eg_i->contig);
				 p_ctg_seq(__func__, eg_j->contig);
				 show_debug_msg(__func__, "Overlapped: %d \n", ol);
				 show_debug_msg(__func__, "Mismatches: %d \n", nm);
				 **/
				// 1. If the overlapping length is shorter than read length,
				// 		We expect that no common reads
				// 2. If the overlapping length is longer than read length,
				//		There mush be some common reads.

				// --------------->
				//           -------------->
				if (ol >= EDGE_OL_THRE && ol * EDGE_OL_PERC > nm) {
					if (ol > hash_o->k || abs(ol - eg_i->len) <= EDGE_OL_THRE
							|| abs(ol - eg_j->len) <= EDGE_OL_THRE) {
						to_merge = 1;
					} else { // Must not share any reads, and there are paired reads
						if (!has_reads_in_common(eg_i, eg_j)) {
							paired_reads = find_unconditional_paired_reads(
									eg_i, eg_j, seqs);
							if (paired_reads->len >= MIN_VALID_PAIRS
									&& check_insert_size(eg_i, eg_j,
											paired_reads, d->insert_size,
											d->sd_insert_size, ol, TAIL_HEAD)) {
								to_merge = 1;
							}
							g_ptr_array_free(paired_reads, TRUE);
							paired_reads = NULL;
						}
					}
					if (to_merge) {
						g_mutex_lock(edge_mutex);
						merge_two_ol_edges(d->rht, d->ht, eg_i, eg_j, ol);
						g_mutex_unlock(edge_mutex);
						eg_i->visited = 0;
						some_one_merged = 1;
						continue;
					} else {
						// ------------->
						//         <--------------
						rev = new_mem_rev_seq(eg_j->contig, eg_j->contig->len,
								0);
						ol = find_ol(eg_i->contig, rev, MAX_EDGE_NM);
						nm = get_mismatches_on_ol(eg_i->contig, rev, ol,
								MAX_EDGE_NM);
						if (ol > hash_o->k || abs(ol - eg_i->len)
								<= EDGE_OL_THRE || abs(ol - eg_j->len)
								<= EDGE_OL_THRE) {
							to_merge = 1;
						} else { // Must not share any reads, and there are paired reads
							if (!has_reads_in_common(eg_i, eg_j)) {
								paired_reads = find_unconditional_paired_reads(
										eg_i, eg_j, seqs);
								if (paired_reads->len >= MIN_VALID_PAIRS
										&& check_insert_size(eg_i, eg_j,
												paired_reads, d->insert_size,
												d->sd_insert_size, ol,
												TAIL_TAIL)) {
									to_merge = 1;
								}
								g_ptr_array_free(paired_reads, TRUE);
								paired_reads = NULL;
							}
						}
						if (to_merge) {
							g_mutex_lock(edge_mutex);
							bwa_free_read_seq(1, eg_j->contig);
							eg_j->contig = rev;
							rev = NULL;
							merge_two_ol_edges(d->rht, d->ht, eg_i, eg_j, ol);
							g_mutex_unlock(edge_mutex);
							eg_i->visited = 0;
							some_one_merged = 1;
							continue; // In case the 'rev' is freed accidently.
						} else {
							// <------------------
							//             -------------->
							bwa_free_read_seq(1, rev);
							rev = new_mem_rev_seq(eg_i->contig,
									eg_i->contig->len, 0);
							ol = find_ol(rev, eg_j->contig, MAX_EDGE_NM);
							nm = get_mismatches_on_ol(rev, eg_j->contig, ol,
									MAX_EDGE_NM);
							if (ol > hash_o->k || abs(ol - eg_i->len)
									<= EDGE_OL_THRE || abs(ol - eg_j->len)
									<= EDGE_OL_THRE) {
								to_merge = 1;
							} else { // Must not share any reads, and there are paired reads
								if (!has_reads_in_common(eg_i, eg_j)) {
									paired_reads
											= find_unconditional_paired_reads(
													eg_i, eg_j, seqs);
									if (paired_reads->len >= MIN_VALID_PAIRS
											&& check_insert_size(eg_i, eg_j,
													paired_reads,
													d->insert_size,
													d->sd_insert_size, ol,
													HEAD_HEAD)) {
										to_merge = 1;
									}
									g_ptr_array_free(paired_reads, TRUE);
									paired_reads = NULL;
								}
							}
							if (to_merge) {
								g_mutex_lock(edge_mutex);
								bwa_free_read_seq(1, eg_i->contig);
								eg_i->contig = rev;
								rev = NULL;
								merge_two_ol_edges(d->rht, d->ht, eg_i, eg_j,
										ol);
								g_mutex_unlock(edge_mutex);
								eg_i->visited = 0;
								some_one_merged = 1;
								continue; // In case the 'rev' is freed accidently.
							}
						} // End of Head-head
						bwa_free_read_seq(1, rev);
						rev = NULL;
					} // End of Tail-tail
				} // End of Tail-head
			} // End of edge candidates
			eg_i->visited = 1;
		}
	}
	if (edge_candidates) {
		g_ptr_array_free(edge_candidates, TRUE);
		edge_candidates = NULL;
	}
	return NULL;
}

void merge_ol_edges(edgearray *single_edges, const int insert_size,
		const int sd_insert_size, hash_table *ht, const int n_threads) {
	int n_per_threads = 0, i = 0;
	merging_paras_t *data;
	GThread *threads[n_threads];
	edge *eg_i = NULL;
	reads_ht *rht = NULL;

	n_per_threads = single_edges->len / n_threads;
	if (!edge_mutex)
		edge_mutex = g_mutex_new();
	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		eg_i->tid = 0;
	}
	data = (merging_paras_t*) calloc(n_threads, sizeof(merging_paras_t));
	show_msg(__func__, "Building hash table for templates...\n");
	rht = build_edges_ht(MATE_OVERLAP_THRE, single_edges);
	show_msg(__func__, "Merging %d templates...\n", single_edges->len);
	for (i = 0; i < n_threads; ++i) {
		data[i].single_edges = single_edges;
		data[i].ht = ht;
		data[i].insert_size = insert_size;
		data[i].start = i * n_per_threads;
		data[i].end = (i + 1) * n_per_threads;
		data[i].tid = i + 1;
		data[i].rht = rht;
		data[i].sd_insert_size = sd_insert_size;
		if (i == n_threads - 1)
			data[i].end = single_edges->len;
		threads[i] = g_thread_create((GThreadFunc) merge_ol_edges_thread,
				data + i, TRUE, NULL);
	}
	for (i = 0; i < n_threads; ++i) {
		g_thread_join(threads[i]);
	}

	for (i = 0; i < single_edges->len; i++) {
		eg_i = g_ptr_array_index(single_edges, i);
		eg_i->visited = 0;
		g_ptr_array_sort(eg_i->reads, (GCompareFunc) cmp_read_by_name);
		if (!eg_i->alive) {
			destroy_eg(eg_i);
			g_ptr_array_remove_index_fast(single_edges, i);
			i--;
		}
	}
	show_msg(__func__, "Merged to %d templates.\n", single_edges->len);
	destroy_reads_ht(rht);
	free(data);
}