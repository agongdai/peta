/*
 * roadmap.c
 *
 *  Created on: Sep 2, 2011
 *      Author: carl
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <glib.h>
#include "roadmap.h"
#include "bwase.h"
#include "peseq.h"
#include "pool.h"
#include "pealn.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"
#include "pechar.h"
#include "utils.h"

int min_edge_len = 0; // Minimum contig length, now is user-specified overlapping len.
// Its index is edge id, say edge A, the value is another edge array, containing edges whose right_ctg is A.
edgearray *egs_right_ctg = 0;

void upd_right_ctg(edge *eg, edge *updated_to, const int change_shift,
		const int add_shift) {
	int i = 0;
	edgearray *right_ctg_connector, *updated_to_right_ctg;
	edge *connector, *right_eg;
	if (!eg)
		return;
	right_ctg_connector = (edgearray*) g_ptr_array_index(egs_right_ctg, eg->id);
	// If updated_to NULL, meaning to remove, so update the egs_right_ctg information
	if (!updated_to) {
		right_eg = eg->right_ctg;
		if (right_eg) {
			right_ctg_connector = (edgearray*) g_ptr_array_index(egs_right_ctg,
					right_eg->id);
			g_ptr_array_remove(right_ctg_connector, eg);
			g_ptr_array_remove(right_eg->in_egs, eg);
		}
		return;
	}
	updated_to_right_ctg = (edgearray*) g_ptr_array_index(egs_right_ctg,
			updated_to->id);
	for (i = 0; i < right_ctg_connector->len; i++) {
		connector = g_ptr_array_index(right_ctg_connector, i);
		if (connector) {
			// Sometimes the input edges eg and updated_to are the same
			// If not the same,update the egs_right_ctg info.
			if (right_ctg_connector != updated_to_right_ctg)
				g_ptr_array_add(updated_to_right_ctg, connector);
			connector->right_ctg = updated_to;
			if (change_shift) {
				// show_debug_msg(__func__, "Shift of edge %d => %d changed: %d => %d \n", connector->id, connector->right_ctg->id, connector->shift, connector->shift + add_shift);
				connector->r_shift += add_shift;
			}
		}
	}
}

/**
 * Contig of left_eg is to be prepend to right_eg
 */
int merge_eg(edge *eg) {
	edge *right_eg, *left_most, *left_eg;
	edgearray *left_in_egs, *out_egs;
	int i = 0;
	if (!eg)
		return 0;
	// For the situation where the left_eg is ori 1, and right_eg is ori 0.
	if (eg->ori && eg->out_egs->len == 1) {
		right_eg = g_ptr_array_index(eg->out_egs, 0);
		if (right_eg->alive && !right_eg->ori && right_eg->in_egs->len == 1) {
			printf("[merge_i] Merging contig %d to %d \n", eg->id,
					right_eg->id);
			merge_seq_to_right(eg->contig, right_eg->contig, 0);
			right_eg->len = right_eg->contig->len;
			right_eg->in_egs = eg->in_egs;
			left_in_egs = eg->in_egs;
			for (i = 0; i < left_in_egs->len; i++) {
				left_most = g_ptr_array_index(left_in_egs, i);
				out_egs = left_most->out_egs;
				g_ptr_array_replace_ptr(out_egs, right_eg, eg);
			}
			eg->alive = 0;
			upd_right_ctg(right_eg, right_eg, 1, eg->len);
			upd_right_ctg(eg, right_eg, 0, 0);
			// Move all reads used by eg to right_eg
			combine_reads(eg, right_eg, 1, 0, 1);
			return 1;
		}
	}
	if (eg->ori) {
		if (eg->in_egs->len != 1)
			return 0;
		left_eg = g_ptr_array_index(eg->in_egs, 0);
		if (left_eg->right_ctg || !eg->alive || left_eg->out_egs->len != 1
				|| !left_eg->alive)
			return 0;
		printf("[merge_i] Merging contig %d to %d \n", left_eg->id, eg->id);
		merge_seq_to_right(left_eg->contig, eg->contig, 0);
		eg->len = eg->contig->len;
		eg->in_egs = left_eg->in_egs;
		left_in_egs = left_eg->in_egs;
		for (i = 0; i < left_in_egs->len; i++) {
			left_most = g_ptr_array_index(left_in_egs, i);
			out_egs = left_most->out_egs;
			g_ptr_array_replace_ptr(out_egs, eg, left_eg);
		}
		left_eg->alive = 0;
		upd_right_ctg(eg, eg, 1, left_eg->len);
		upd_right_ctg(left_eg, eg, 0, 0);
		combine_reads(left_eg, eg, 1, 0, 1);
		return 1;
	} else {
		if (eg->out_egs->len != 1)
			return 0;
		right_eg = g_ptr_array_index(eg->out_egs, 0);
		if (eg->right_ctg || !eg->alive || right_eg->in_egs->len != 1
				|| !right_eg->alive)
			return 0;
		printf("[merge_i] Merging contig %d to %d \n", eg->id, right_eg->id);
		merge_seq_to_right(eg->contig, right_eg->contig, 0);
		right_eg->len = right_eg->contig->len;
		right_eg->in_egs = eg->in_egs;
		left_in_egs = eg->in_egs;
		for (i = 0; i < left_in_egs->len; i++) {
			left_most = g_ptr_array_index(left_in_egs, i);
			out_egs = left_most->out_egs;
			g_ptr_array_replace_ptr(out_egs, right_eg, eg);
		}
		eg->alive = 0;
		upd_right_ctg(right_eg, right_eg, 1, eg->len);
		upd_right_ctg(eg, right_eg, 0, 0);
		combine_reads(eg, right_eg, 1, 0, 1);
		return 1;
	}
}

/**
 * left_in_egs         left_eg         out_egs
 * -----------                         -----------
 * -----------                         -----------
 * -----------         =========       -----------
 * -----------                         -----------
 * -----------
 *
 * left_in_egs            left_eg   out_egs
 * -----------           =========-----------
 * -----------           =========-----------
 * -----------           =========-----------
 * -----------           =========-----------
 * -----------
 */
int backward_branches(edge *left_eg) {
	bwa_seq_t *left_ctg;
	edge *eg_i = 0, *left_most;
	edgearray *left_in_egs, *out_egs, *egs;
	int i = 0, j = 0;
	if (!left_eg || left_eg->right_ctg || left_eg->len >= min_edge_len
			|| !left_eg->alive)
		return 0;
	left_in_egs = left_eg->in_egs;
	out_egs = left_eg->out_egs;
	//
	if (out_egs->len <= 1)
		return 0;
	left_ctg = left_eg->contig;

	egs = (edgearray*) g_ptr_array_index(egs_right_ctg, left_eg->id);
	// Will not backward if some edge's right_ctg is left_eg, because after backwarding,
	//    cannot update the right_ctg info.
	if (egs->len > 0)
		return 0;

	printf("[backward_branches] Backwarding branches of edge %d... \n",
			left_eg->id);
	show_debug_msg(__func__, "# of out_egs: %d\n", out_egs->len);
	for (i = 0; i < out_egs->len; i++) {
		eg_i = g_ptr_array_index(out_egs, i);
//		show_debug_msg(__func__, "Pointer of out_egs: %p\n", out_egs);
//		show_debug_msg(__func__, "# of out_egs: %d\n", out_egs->len);
//		show_debug_msg(__func__, "Out edges of %d: %d\n", left_eg->id, i);
//		p_flat_eg(eg_i);
		printf("\n");
		if (!eg_i->contig->seq)
			continue;
		merge_seq_to_right(left_ctg, eg_i->contig, 0);
		eg_i->len = eg_i->contig->len;
		eg_i->in_egs = left_eg->in_egs;
		for (j = 0; j < left_in_egs->len; j++) {
			left_most = g_ptr_array_index(left_in_egs, j);
			egs = left_most->out_egs;
			// Out edges of eg_i should be appended to left_most correctly
			// Cannot be just left_most->out_egs = out_egs, because it may contain some other edges.
			if (i == 0)
				g_ptr_array_replace_ptr(egs, eg_i, left_eg);
			else
				g_ptr_array_uni_add(egs, eg_i);
		}
		// Move all used reads of left_eg to eg_i, updating the shift value
		combine_reads(left_eg, eg_i, 1, 0, 1);
	}
	upd_right_ctg(eg_i, eg_i, 1, left_eg->len);
	upd_right_ctg(left_eg, eg_i, 0, 0);
	left_eg->alive = 0;
	return 1;
}

/**
 * Return value indicates whether collapsed (1) or not (0)
 */
int collapse_branches(edge *left_eg) {
	edgearray *out_egs, *distinct_egs, *egs, *collapsed_right;
	int i = 0, j = 0, exist_index = NOT_FOUND;
	edge *eg, *collapsed_to = 0, *tmp;
	if (!left_eg || left_eg->right_ctg)
		return 0;
	out_egs = left_eg->out_egs;
	if (out_egs->len < 2)
		return 0;
	distinct_egs = g_ptr_array_sized_new(out_egs->len);
	// distinct_edges stores all distinct edges of the out edges
	for (i = 0; i < out_egs->len; i++) {
		eg = g_ptr_array_index(out_egs, i);
		if (eg->len >= min_edge_len * 2) {
			g_ptr_array_add(distinct_egs, eg);
		} else {
			exist_index = edgearray_find_similar(distinct_egs, eg);
			if (exist_index == NOT_FOUND)
				g_ptr_array_add(distinct_egs, eg);
			else
				// !!! Temp storage, indicating the edge to be merged to
				eg->right_ctg = g_ptr_array_index(distinct_egs, exist_index);
		}
	}
	// TODO: local alignment to determine what to retain,
	//		currently just longest survives
	// If there are similar out edges
	if (distinct_egs->len < out_egs->len) {
		printf("[collapse_branches] Collapsing branches of edge %d... \n",
				left_eg->id);
		for (i = 0; i < out_egs->len; i++) {
			eg = g_ptr_array_index(out_egs, i);
			exist_index = edgearray_find(distinct_egs, eg);
			if (exist_index != NOT_FOUND || !eg->right_ctg)
				continue;

			egs = eg->out_egs;
			collapsed_to = eg->right_ctg; // Now the value is just temp!!!
			assert(collapsed_to && collapsed_to->id >= 0);
			show_debug_msg(__func__,
					"Collapsing branches of edge %d: edge %d => %d... \n",
					left_eg->id, eg->id, collapsed_to->id);
			// If right_ctg value is not NULL, collapsed_to and eg must have the same value
			// Here remove the edge to be collapsed.
			if (collapsed_to->right_ctg) {
				g_ptr_array_remove(collapsed_to->right_ctg->in_egs, eg);
			}
			collapsed_right = collapsed_to->out_egs;
			// The two edges may connect to the same node
			if (egs != collapsed_right) {
				// For each edge to be collapsed, repointing its out edges to the similar distinct edge.
				for (j = 0; j < egs->len; j++) {
					tmp = g_ptr_array_index(egs, j);
					if (edgearray_find(collapsed_right, tmp) == NOT_FOUND)
						g_ptr_array_add(collapsed_right, tmp);
					g_ptr_array_replace_ptr(tmp->in_egs, collapsed_to, eg);
				}
			} else {
				// eg's right_ctg is collapsed_to, so remove from the in-edges of their out-edegs.
				for (j = 0; j < egs->len; j++) {
					tmp = g_ptr_array_index(egs, j);
					g_ptr_array_remove(tmp->in_egs, eg);
				}
			}
			upd_right_ctg(eg, collapsed_to, 1, 0);
			eg->alive = 0;
			// Add all used reads of eg to collapsed_to,
			// Not changing the shift value
			combine_reads(collapsed_to, eg, 0, 0, 0);
		}

		// Clear all edges, and readd all distinct edges
		g_ptr_array_remove_range(out_egs, 0, out_egs->len);
		g_ptr_array_concat(out_egs, distinct_egs);
		g_ptr_array_free(distinct_egs, TRUE);
		return 1;
	} else
		return 0;
}

roadmap *new_rm() {
	roadmap *rm = (roadmap*) malloc(sizeof(roadmap));
	rm->n_node = rm->start_eg_n = 0;
	rm->pet_e = rm->pet_s = 0;
	rm->start_egs = g_ptr_array_sized_new(INIT_PET_N);
	return rm;
}

edge *new_eg() {
	edge *eg = (edge*) malloc(sizeof(edge));
	eg->contig = 0;
	eg->in_egs = g_ptr_array_sized_new(MAX_N_EDGE_IN);
	eg->out_egs = g_ptr_array_sized_new(MAX_N_EDGE_OUT);
	eg->reads = g_ptr_array_sized_new(INIT_N_READ_USED);
	eg->name = NULL;
	eg->right_ctg = NULL;
	eg->left_ctg = NULL;
	eg->len = 0;
	eg->r_shift = 0;
	eg->l_shift = 0;
	eg->id = 0;
	eg->visited = 0;
	eg->alive = 1;
	eg->is_root = 0;
	eg->ori = 0;
	eg->gaps = g_ptr_array_sized_new(INIT_N_GAPS);
	return eg;
}

void free_readarray(readarray *ra) {
	if (!ra)
		return;
	g_ptr_array_free(ra, TRUE);
}

void destroy_eg(edge *eg) {
	eg_gap *gap = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	if (eg) {
		bwa_free_read_seq(1, eg->contig); // bug if free it
		g_ptr_array_free(eg->in_egs, TRUE);
		if (!eg->right_ctg) {
			// If eg's right contig is not null, its out_egs is set to be right contig's out_egs
			g_ptr_array_free(eg->out_egs, TRUE);
		}
		for (i = 0; i < eg->reads->len; i++) {
			read = g_ptr_array_index(eg->reads, i);
			read->used = 0;
			read->contig_id = UNUSED_CONTIG_ID;
		}
		while (eg->reads->len > 0) {
			g_ptr_array_remove_index_fast(eg->reads, 0);
		}
		free_readarray(eg->reads);
		for (i = 0; i < eg->gaps->len; i++) {
			gap = g_ptr_array_index(eg->gaps, i);
			free_eg_gap(gap);
		}
		free_readarray(eg->gaps);
		eg->alive = 0;
		free(eg->name);
		free(eg);
	}
}

/**
 * The reads used by this edge could be reused later (set attribute 'used' to be false).
 */
void free_eg(edge *eg, const int ori) {
	eg_gap *gap = NULL;
	int i = 0;
	bwa_seq_t *read = NULL;
	show_debug_msg("REMOVED", "Edge [%d, %d] removed \n", eg->id, eg->len);
	if (eg) {
		bwa_free_read_seq(1, eg->contig);
		if (eg->right_ctg) {
			g_ptr_array_remove(eg->right_ctg->in_egs, eg);
		}
		if (ori)
			free_readarray(eg->in_egs);
		if (!eg->right_ctg && !ori) {
			// If eg's right contig is not null, its out_egs is set to be right contig's out_egs
			free_readarray(eg->out_egs);
		}
		eg->out_egs = NULL;
		for (i = 0; i < eg->reads->len; i++) {
			read = g_ptr_array_index(eg->reads, i);
			read->used = 0;
			read->contig_id = UNUSED_CONTIG_ID;
		}
		while (eg->reads->len > 0) {
			g_ptr_array_remove_index_fast(eg->reads, 0);
		}
		free_readarray(eg->reads);
		for (i = 0; i < eg->gaps->len; i++) {
			gap = g_ptr_array_index(eg->gaps, i);
			free_eg_gap(gap);
		}
		free_readarray(eg->gaps);
		eg->alive = 0;
		eg->right_ctg = NULL;
		free(eg->name);
		free(eg);
	}
}

void cut_connection(edge *ass_eg, edge *tmp_eg, const int ori) {
	if (ori) {
		g_ptr_array_remove(ass_eg->in_egs, tmp_eg);
	} else {
		g_ptr_array_remove(ass_eg->out_egs, tmp_eg);
	}
}

void free_branch(edge *eg, const int ori, edgearray *all_edges) {
	edgearray *children = NULL;
	int i = 0;
	edge *child = NULL;
	if (!eg)
		return;
	g_ptr_array_remove(all_edges, eg);
	children = ori ? eg->in_egs : eg->out_egs;
	if (!ori && eg->right_ctg) {
	} else {
		for (i = 0; i < children->len; i++) {
			child = g_ptr_array_index(children, i);
			free_branch(child, ori, all_edges);
		}
	}
	free_eg(eg, ori);
}

/**
 * in_egs              left_eg             out_egs
 * -----------                                 -----------
 * -----------                                 -----------
 * -----------         =========   -----------
 * -----------                                 -----------
 * -----------
 *
 * in_egs                out_egs
 * -----------           -----------
 * -----------           -----------
 * -----------           -----------
 * -----------           -----------
 * -----------
 */
int rm_eg(edge *eg) {
	int i = 0, j = 0;
	edge *eg_i, *tmp;
	edgearray *out_egs, *in_egs, *egs;
	if (!eg || !eg->alive)
		return 0;
	printf("[rm_eg] Removing Contig %d... \n", eg->id);
	out_egs = eg->out_egs;
	in_egs = eg->in_egs;
	// For each out-edge, remove current edge from its in-edge list.
	for (i = 0; i < out_egs->len; i++) {
		eg_i = g_ptr_array_index(out_egs, i);
		egs = eg_i->in_egs;
		for (j = 0; j < egs->len; j++) {
			tmp = g_ptr_array_index(egs, j);
			if (tmp->id == eg->id) {
				g_ptr_array_remove_index(egs, j);
				break;
			}
		}
	}
	// For each in-edge, remove current edge from its out-edge list.
	for (i = 0; i < in_egs->len; i++) {
		eg_i = g_ptr_array_index(in_egs, i);
		egs = eg_i->out_egs;
		for (j = 0; j < egs->len; j++) {
			tmp = g_ptr_array_index(egs, j);
			if (tmp->id == eg->id) {
				g_ptr_array_remove_index(egs, j);
				break;
			}
		}
	}
	if (eg->right_ctg) {
		in_egs = eg->right_ctg->in_egs;
		for (i = 0; i < in_egs->len; i++) {
			eg_i = g_ptr_array_index(in_egs, i);
			if (eg_i->id == eg->id) {
				g_ptr_array_remove_index(in_egs, i);
				break;
			}
		}
	}
	eg->alive = 0;
	clear_used_reads(eg, 1);
	upd_right_ctg(eg, 0, 0, 0);
	return 1;
}

int prune_eg(edge *eg) {
	edgearray *in_egs, *out_egs, *sbls;
	edge *e, *r_ctg, *eg_i;
	int i = 0;
	if (!eg || !eg->alive)
		return 0;

	in_egs = eg->in_egs;
	out_egs = eg->out_egs;

	if (eg->len < MINCONTIG) {
		if (eg->ori) {
			sbls = get_sbls(eg);
			if (in_egs->len == 0 && sbls->len != 0)
				return rm_eg(eg);
		} else {
			if (!eg->right_ctg) {
				// If is a leaf
				sbls = get_sbls(eg);
				if (out_egs->len == 0 && sbls->len != 0) {
					return rm_eg(eg);
				}
			} else {
				// If the current edge "hangs" on some contig and length not long enough.
				//  |-> eg: ----()---------------
				//  |-> eg->right_ctg (shift = 4)
				if (is_sbl(eg, eg->right_ctg)
						&& ((eg->len - eg->r_shift) < MINCONTIG
								|| (abs(
										(eg->len + eg->r_shift
												- eg->right_ctg->len))
										< MINCONTIG))) {
					return rm_eg(eg);
				}
				//                  |-> eg_0: aaaacccgg |-> e
				// sbls:            |-> eg_1
				//                  eg: aaaaccc (eg->right_ctg = e, shift = -2)
				if (eg->len < MINCONTIG && eg->right_ctg
						&& eg->r_shift < MINCONTIG && in_egs->len > 0) {
					sbls = get_sbls(eg);
					for (i = 0; i < sbls->len; i++) {
						eg_i = g_ptr_array_index(sbls, i);
						if (has_edge(eg_i, eg->right_ctg->id, 0)) {
							if (abs(
									eg_i->len
											- (eg->len - eg->r_shift)) < MINCONTIG) {
								return rm_eg(eg);
							}
						}
					}
				}
			}
			if (out_egs->len == 1) {
				// |-> r_ctg: ----()----------------
				// |
				// |-> eg: ----- |-> e: -------- (e->right_ctg = r_ctg, shift = 4)
				e = g_ptr_array_index(out_egs, 0);
				r_ctg = e->right_ctg;
				if (r_ctg && is_sbl(eg, r_ctg) && (eg->len + e->len
						- e->r_shift) < MINCONTIG) {rm_eg(e);
				return rm_eg(eg);
			}
		}
	}
}
	return 0;
}

int post_pro_edges(edgearray *all_edges) {
	int updated = 0, i = 0;
	edge *eg, *right_ctg;
	edgearray *right_ctg_connectors;
	for (i = all_edges->len - 1; i >= 0; i--) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg && eg->alive) {
			// show_debug_msg(__func__, "Processing edge %d...\n", eg->id);
			right_ctg = eg->right_ctg;
			right_ctg_connectors = g_ptr_array_index(egs_right_ctg, eg->id);
			if (right_ctg_connectors->len == 0)
				updated = updated | prune_eg(eg);
			if (eg->alive)
				updated = updated | merge_eg(eg);
		}
	}
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (eg && eg->alive) {
			// show_debug_msg(__func__, "Processing edge %d...\n", eg->id);
			updated = updated | collapse_branches(eg);
			updated = updated | backward_branches(eg);
		}
	}
	return updated;
}

void post_pro(edgearray *all_edges, const ass_opt *opt) {
	int i = 0, updated;
	edge *eg, *right_ctg;
	edgearray *right_ctg_connectors;
	min_edge_len = opt->ol;
	egs_right_ctg = g_ptr_array_sized_new(all_edges->len + 1);
	right_ctg_connectors = g_ptr_array_sized_new(RIGHT_CONTIG_CONNETER_COUNT);
	char *graph_fn = malloc(BUFSIZE);

	// egs_right_ctg: 1 => 2,3,4, meaning that edge 2, 3, and 4's right_ctg is edge 1
	for (i = 0; i <= all_edges->len; i++) {
		right_ctg_connectors = g_ptr_array_sized_new(
				RIGHT_CONTIG_CONNETER_COUNT);
		g_ptr_array_add(egs_right_ctg, right_ctg_connectors);
	}

	for (i = 0; i < all_edges->len; i++) {
		// Store length of all original contigs.
		eg = g_ptr_array_index(all_edges, i);
		if (eg) {
			right_ctg = eg->right_ctg;
			// Store all edges which is the right_ctg of some other edges.
			if (right_ctg) {
				right_ctg_connectors = (edgearray*) g_ptr_array_index(
						egs_right_ctg, right_ctg->id);
				g_ptr_array_add(right_ctg_connectors, eg);
			}
		}
	}

	// Loop until cannot be updated further
	i = 0;
	do {
		updated = post_pro_edges(all_edges);
		sprintf(graph_fn, "graph/rm_after_update_%d.dot", i++);
		graph_by_edges(all_edges, graph_fn);
	} while (updated);
	free(graph_fn);
}
