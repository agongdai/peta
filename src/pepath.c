/*
 * pepath.c
 *
 *  Created on: Sep 6, 2011
 *      Author: Cai Shaojiang
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <glib.h>
#include "pepath.h"
#include "roadmap.h"
#include "bwase.h"
#include "pechar.h"
#include "utils.h"
#include "pehash.h"
#include "rnaseq.h"
#include "peseq.h"
#include "readrm.h"
#include "edge.h"
#include "edgelist.h"
#include "scaffolding.h"

int path_id = 0;

extern unsigned char nst_nt4_table[256];

rm_path *new_path() {
	rm_path *p = (rm_path*) malloc(sizeof(rm_path));
	p->edges = g_ptr_array_sized_new(16);
	p->n_ctgs = 0;
	p->len = 0;
	p->alive = 1;
	p->seq = NULL;
	p->reads = g_ptr_array_sized_new(16);
	return p;
}

void p_path(const rm_path *p) {
	int i = 0;
	edge *eg = NULL;
	if (!p) {
		printf("Empty path! \n");
		return;
	}
	printf("[p_path] ----------------------------------------\n");
	printf("[p_path] Path %d (%p): %d \n", p->id, p, p->len);
	for (i = 0; i < p->edges->len; i++) {
		eg = g_ptr_array_index(p->edges, i);
		printf("[p_path: ori %d] \t %d: Contig [%d: %d] (%d, %d)\n", eg->ori,
				i, eg->id, eg->len, eg->right_ctg ? eg->right_ctg->id : 0,
				eg->r_shift);
	}
	p_ctg_seq("Path seq", p->seq);
	printf("[p_path] ----------------------------------------\n");
}

/**
 * Destroy a roadmap path
 */
void destroy_path(rm_path *p) {
	if (p) {
		g_ptr_array_free(p->reads, TRUE);
		g_ptr_array_free(p->edges, TRUE);
		bwa_free_read_seq(1, p->seq);
		free(p);
	}
}

void save_paths(GPtrArray *paths, const char *tx_fn, const int min_len) {
	char header[BUFSIZE];
	int i = 0;
	FILE *tx = xopen(tx_fn, "w");
	rm_path *p = NULL;
	show_msg(__func__, "Saving transcripts to %s...\n", tx_fn);
	for (i = 0; i < paths->len; i++) {
		p = g_ptr_array_index(paths, i);
		//p_path(p);
		if (p->len >= min_len) {
			sprintf(header, ">%d len=%d \n", i, p->len);
			save_con(header, p->seq, tx);
		}
	}
	fclose(tx);
}

void sync_path(rm_path *p) {
	bwa_seq_t *seq = NULL, *rev_seq = NULL;
	int i = 0, shift = 0;
	edge *eg = NULL, *eg_next = NULL;
	p->n_ctgs = p->edges->len;
	seq = blank_seq();
	if (p->n_ctgs > 0) {
		for (i = 0; i < p->n_ctgs; i++) {
			eg = g_ptr_array_index(p->edges, i);
			if (shift < eg->len) {
				if (shift < 0)
					shift = 0;
				if (eg->ori == 1) {
					rev_seq = new_mem_rev_seq(eg->contig, eg->len, 0);
					merge_seq(seq, rev_seq, shift);
				} else {
					merge_seq(seq, eg->contig, shift);
				}
			} else {
				if (i + 1 < p->n_ctgs) {
					eg_next = g_ptr_array_index(p->edges, i + 1);
					if (eg_next->r_shift == 0)
						eg->r_shift = shift - eg->len;
				}
			}
			shift = eg->r_shift;
		}
	}
	bwa_free_read_seq(1, p->seq);
	p->seq = seq;
	p->len = seq->len;
}

/**
 * Add an edge to a path; if adding to head, create a brand new pointer array
 */
void add_edge_to_path(rm_path *path, edge *eg, const int head) {
	GPtrArray *edges = NULL;
	int i = 0, len = 0;
	edge *eg_i = NULL;
	if (head) {
		edges = g_ptr_array_sized_new(path->n_ctgs + 4);
		g_ptr_array_add(edges, eg); // Put the new one as the first one.
		len += eg->len;
		for (i = 0; i < path->edges->len; i++) {
			eg_i = g_ptr_array_index(path->edges, i);
			g_ptr_array_add(edges, eg_i);
			len += eg_i->len;
		}
		path->n_ctgs = edges->len;
		path->len = len;
		g_ptr_array_free(path->edges, TRUE); // Free the old array
		path->edges = edges;
	} else {
		g_ptr_array_add(path->edges, eg);
		path->len += eg->len;
		path->n_ctgs++;
	}
}

rm_path *get_single_edge_path(edge *eg) {
	rm_path *p = new_path();
	p->id = eg->id;
	add_edge_to_path(p, eg, 0);
	return p;
}

/**
 * Clone a path
 */
rm_path *clone_path(const rm_path *p) {
	int i = 0;
	edge *eg = NULL;
	rm_path *new_p = NULL;

	new_p = new_path();
	new_p->edges = g_ptr_array_sized_new(p->n_ctgs);
	if (!p)
		return new_p;
	for (i = 0; i < p->n_ctgs; i++) {
		eg = g_ptr_array_index(p->edges, i);
		g_ptr_array_add(new_p->edges, eg);
	}
	new_p->n_ctgs = new_p->edges->len;
	new_p->len = p->len;
	new_p->alive = 1;
	new_p->id = p->id;
	return new_p;
}

/**
 * Fork a path from the middle
 */
rm_path *fork_path(const rm_path *p, edge *right_connect, const int joint_point) {
	rm_path *forked = NULL;
	edge *eg_i = NULL;
	int i = 0;
	forked = new_path();
	add_edge_to_path(forked, right_connect, 0);
	if (joint_point >= 0 && joint_point < p->n_ctgs) {
		for (i = joint_point; i < p->n_ctgs; i++) {
			eg_i = g_ptr_array_index(p->edges, i);
			add_edge_to_path(forked, eg_i, 0);
		}
	}
	forked->id = p->id;
	return forked;
}

int same_paths(rm_path *p1, rm_path *p2) {
	int i = 0;
	edge *eg_1 = NULL, *eg_2 = NULL;
	if (!p1 || !p2 || p1->n_ctgs != p2->n_ctgs)
		return 0;
	for (i = 0; i < p1->n_ctgs; i++) {
		eg_1 = g_ptr_array_index(p1->edges, i);
		eg_2 = g_ptr_array_index(p2->edges, i);
		if (eg_1 != eg_2)
			return 0;
	}
	return 1;
}

/**
 * Determine whether a path is a subpath of another.
 * For example: p1 = [e1, e2, e3], p2 = [e2, e3]
 * Concatenate edge ids: p1 = 'e1,e2,e3,', p2 = 'e2,e3,'
 * p2 string is as substring of p1, so return 1.
 */
int is_sub_path(rm_path *p1, rm_path *p2) {
	int id_str_len = 0, i = 0;
	edge *eg_i = 0;
	char str_id[32], *ids_1 = NULL, *ids_2 = NULL, *is_sub = NULL;
	if (p1->n_ctgs < p2->n_ctgs)
		return 0;
	for (i = 0; i < p1->n_ctgs; i++) {
		eg_i = g_ptr_array_index(p1->edges, i);
		sprintf(str_id, "%d,", eg_i->id);
		id_str_len += strlen(str_id);
	}
	ids_1 = (char*) calloc(id_str_len + 1, sizeof(char));
	for (i = 0; i < p1->n_ctgs; i++) {
		eg_i = g_ptr_array_index(p1->edges, i);
		sprintf(str_id, "%d,", eg_i->id);
		strcat(ids_1, str_id);
	}
	for (i = 0; i < p2->n_ctgs; i++) {
		eg_i = g_ptr_array_index(p2->edges, i);
		sprintf(str_id, "%d,", eg_i->id);
		id_str_len += strlen(str_id);
	}
	ids_2 = (char*) calloc(id_str_len + 1, sizeof(char));
	for (i = 0; i < p2->n_ctgs; i++) {
		eg_i = g_ptr_array_index(p2->edges, i);
		sprintf(str_id, "%d,", eg_i->id);
		strcat(ids_2, str_id);
	}
	is_sub = strstr(ids_1, ids_2);
	free(ids_1);
	free(ids_2);
	if (is_sub == NULL)
		return 0;
	else
		return 1;
}

int find_path(GPtrArray *paths, rm_path *p) {
	int i = 0;
	rm_path *path = NULL;
	for (i = 0; i < paths->len; i++) {
		path = g_ptr_array_index(paths, i);
		if (same_paths(path, p))
			return i;
	}
	return -1;
}

void add_uni_path(GPtrArray *paths, rm_path *p) {
	if (find_path(paths, p) == -1) {
		g_ptr_array_add(paths, p);
	}
}

/*
 * Recursive method to get all edges reachable by the input edge
 *
 * Block level indicates how many levels in the graph
 */
void get_block_edges(edge *eg, GPtrArray *block) {
	edge *in_out_eg = NULL;
	int i = 0;
	if (!eg->visited) {
		g_ptr_array_add(block, eg);
		eg->visited = 1;
	}
	for (i = 0; i < eg->in_egs->len; i++) {
		in_out_eg = g_ptr_array_index(eg->in_egs, i);
		if (in_out_eg->alive && !in_out_eg->visited) {
			g_ptr_array_add(block, in_out_eg);
			in_out_eg->visited = 1;
			in_out_eg->level = eg->level - 1;
			get_block_edges(in_out_eg, block);
		}
	}
	for (i = 0; i < eg->out_egs->len; i++) {
		in_out_eg = g_ptr_array_index(eg->out_egs, i);
		if (in_out_eg->alive && !in_out_eg->visited) {
			g_ptr_array_add(block, in_out_eg);
			in_out_eg->visited = 1;
			in_out_eg->level = eg->level + 1;
			get_block_edges(in_out_eg, block);
		}
	}
	if (eg->right_ctg && eg->right_ctg->alive && !eg->right_ctg->visited) {
		g_ptr_array_add(block, eg->right_ctg);
		eg->right_ctg->visited = 1;
		get_block_edges(eg->right_ctg, block);
	}
}

/**
 * For a list of paths, try to attach an edge to them.
 *
 * Only if the edge could be added to the head or tail of a path
 */
GPtrArray *attach_edge_to_paths(GPtrArray *paths, edge *eg) {
	int i = 0, j = 0;
	rm_path *p = NULL, *new_p = NULL;
	edge *head = NULL, *tail = NULL, *eg_i = NULL;
	edgearray *edges = NULL;
	GPtrArray *paths_to_add = NULL;

	paths_to_add = g_ptr_array_sized_new(8);
	if (!paths || paths->len == 0 || !eg)
		return paths_to_add;
	for (i = 0; i < paths->len; i++) {
		p = g_ptr_array_index(paths, i);
		edges = p->edges;
		head = g_ptr_array_index(edges, 0);
		tail = g_ptr_array_index(edges, edges->len - 1);
		if (edgearray_find(tail->out_egs, eg) != -1) {
			new_p = clone_path(p);
			add_edge_to_path(new_p, eg, 0);
			p->alive = 0;
			add_uni_path(paths_to_add, new_p);
		} else {
			if (edgearray_find(head->in_egs, eg) != -1) {
				new_p = clone_path(p);
				add_edge_to_path(new_p, eg, 1);
				p->alive = 0;
				add_uni_path(paths_to_add, new_p);
			} else {
				for (j = 0; j < p->n_ctgs; j++) {
					eg_i = g_ptr_array_index(p->edges, j);
					if (eg->right_ctg == eg_i) {
						new_p = fork_path(p, eg, j);
						add_uni_path(paths_to_add, new_p);
					}
				}
			}
		}
	}

	return paths_to_add;
}

int has_super_path(GPtrArray *paths, rm_path *p) {
	int i = 0;
	rm_path *p_i = NULL;
	for (i = 0; i < paths->len; i++) {
		p_i = g_ptr_array_index(paths, i);
		if (p_i != p && is_sub_path(p_i, p))
			return 1;
	}
	return 0;
}

/*
 * Get all paths in the block
 */
GPtrArray *iterate_block(GPtrArray *block) {
	edge *eg = NULL;
	rm_path *p = NULL;
	int i = 0, has_fresh = 1;
	edgearray *edges_in_paths = NULL;
	GPtrArray *paths = NULL, *block_paths_to_add = NULL;

	paths = g_ptr_array_sized_new(8);
	// Add the root edges, each of which is a stand-alone path
	if (block->len == 1) {
		eg = g_ptr_array_index(block, 0);
		p = get_single_edge_path(eg);
		g_ptr_array_add(paths, p);
		return paths;
	}
	edges_in_paths = g_ptr_array_sized_new(block->len + 1);
	for (i = 0; i < block->len; i++) {
		eg = g_ptr_array_index(block, i);
		if (eg->is_root) {
			p = get_single_edge_path(eg);
			g_ptr_array_add(paths, p);
			g_ptr_array_add(edges_in_paths, eg);
		}
	}

	// If there is no edge which is able to be added, stop
	while (has_fresh) {
		has_fresh = 0;
		for (i = 0; i < block->len; i++) {
			eg = g_ptr_array_index(block, i);
			if (edgearray_find(edges_in_paths, eg) == NOT_FOUND) {
				block_paths_to_add = attach_edge_to_paths(paths, eg); // All paths in the block
				if (block_paths_to_add->len > 0) {
					has_fresh = 1;
					g_ptr_array_concat(paths, block_paths_to_add); // add new paths
					g_ptr_array_add(edges_in_paths, eg); // Mark this edge has been visited
				}
				g_ptr_array_free(block_paths_to_add, TRUE);
			}
		}
	}

	// Remove the paths which are not alive already.
	for (i = 0; i < paths->len; i++) {
		p = g_ptr_array_index(paths, i);
		// If a path is a sub-path of another, remove it
		if (!p->alive || has_super_path(paths, p)) {
			g_ptr_array_remove_index_fast(paths, i);
			i--;
			destroy_path(p);
		}
	}

	// For cases: eg->right_ctg = right_eg, with shift value 500.
	// Paths could be: eg, right_eg[500:-1].
	// In this case, the starting 500 bases of right_eg is not considered at all.
	// So we add these 'right_eg's as a single-edge path
	g_ptr_array_free(edges_in_paths, TRUE);
	edges_in_paths = g_ptr_array_sized_new(block->len + 1);
	for (i = 0; i < block->len; i++) {
		eg = g_ptr_array_index(block, i);
		if (eg->right_ctg && edgearray_find(edges_in_paths, eg->right_ctg)
				== NOT_FOUND) {
			p = get_single_edge_path(eg->right_ctg);
			g_ptr_array_add(paths, p);
			g_ptr_array_add(edges_in_paths, eg->right_ctg);
		}
	}

	g_ptr_array_free(edges_in_paths, TRUE);
	return paths;
}

/**
 * Get the max difference of levels in the block paths
 */
int get_level_n(GPtrArray *block) {
	int min_level = 0, max_level = 0;
	int i = 0;
	edge *eg = NULL;
	for (i = 0; i < block->len; i++) {
		eg = g_ptr_array_index(block, i);
		if (eg->level < min_level)
			min_level = eg->level;
		if (eg->level > max_level)
			max_level = eg->level;
	}
	return max_level - min_level;
}

/**
 * Check whether there is counter level
 */
int has_possible_circles(GPtrArray *block) {
	int i = 0;
	edge *eg = NULL, *eg_right = NULL;
	for (i = 0; i < block->len; i++) {
		eg = g_ptr_array_index(block, i);
		if (eg->right_ctg) {
			eg_right = eg->right_ctg;
			if (eg_right->level <= eg->level)
				return 1;
		}
	}
	return 0;
}

GPtrArray *get_single_block_paths(GPtrArray *block) {
	GPtrArray *single_paths = NULL;
	int i = 0;
	rm_path *p;
	edge *eg = NULL;

	single_paths = g_ptr_array_sized_new(block->len);
	for (i = 0; i < block->len; i++) {
		eg = g_ptr_array_index(block, i);
		if (eg->len >= MIN_TX_LEN) {
			p = new_path();
			add_edge_to_path(p, eg, 0);
			p->id = eg->id;
			g_ptr_array_add(single_paths, p);
		}
	}
	return single_paths;
}

/**
 * Mark the duplicate edges
 */
void mark_duplicate_edges(edgearray *block) {
	edge *eg_i = NULL, *eg_j = NULL;
	int i = 0, j = 0;
	for (i = 0; i < block->len; i++) {
		eg_i = g_ptr_array_index(block, i);
		for (j = 0; j < block->len; j++) {
			eg_j = g_ptr_array_index(block, j);
			if (eg_i != eg_j && eg_j->alive && similar_seqs(eg_i->contig,
					eg_j->contig, PATH_MISMATCHES, PATH_MAX_GAPS, SCORE_MATCH,
					SCORE_MISMATCH, SCORE_GAP)) {
				if (eg_i->right_ctg == eg_j->right_ctg && abs(eg_i->r_shift
						- eg_j->r_shift) <= PATH_MISMATCHES) {
					eg_i->alive = 0;
					show_debug_msg("DUPLICATE", "[%d: %d] <=> [%d, %d] \n",
							eg_i->id, eg_i->len, eg_j->id, eg_j->len);
					break;
				}
			}
		}
	}
}

void get_path_ori(rm_path *path, bwa_seq_t *seqs) {
	edge *eg_left = NULL, *eg_right = NULL;
	int i = 0;
	if (path->edges->len > 1) {
		for (i = 0; i < path->edges->len - 1; i++) {
			eg_left = g_ptr_array_index(path->edges, i);
			eg_right = g_ptr_array_index(path->edges, i + 1);
			get_edges_ori(eg_left, eg_right, seqs);
		}
	}
}

void determine_paths_ori(GPtrArray *paths, bwa_seq_t *seqs) {
	int i = 0;
	rm_path *p = NULL;
	for (i = 0; i < paths->len; i++) {
		p = g_ptr_array_index(paths, i);
		get_path_ori(p, seqs);
		sync_path(p);
	}
}

void mark_duplicate_paths(GPtrArray *paths) {
	int i = 0, j = 0, similarity_score = 0, report_unit = 0;
	rm_path *path_i = NULL, *path_j = NULL;
	edge *eg = NULL;
	show_debug_msg(__func__, "Removing paths which have duplicate edges...\n");
	report_unit = paths->len / 10;
	// Mark a path as not alive if some edge on it is not alive
	for (i = 0; i < paths->len; i++) {
		if (report_unit >= 20 && (i) % report_unit == 0)
			show_msg(__func__, "Progress 1/2: %d/%d...\n", i, paths->len);
		path_i = g_ptr_array_index(paths, i);
		// If any edge in the path is not alive, remove this path
		for (j = 0; j < path_i->n_ctgs; j++) {
			eg = g_ptr_array_index(path_i->edges, j);
			if (!eg->alive) {
				path_i->alive = 0;
				break;
			}
		}
		if (!path_i->alive) {
			g_ptr_array_remove_fast(paths, path_i);
			i--;
		} else {
			sync_path(path_i); // Populate the path sequence
		}
	}
	show_debug_msg(__func__,
			"Smith-waterman algorithm to remove duplicate paths...\n");
	// Mark a path as not alive if some alive path is similar to it.
	for (i = 0; i < paths->len; i++) {
		if (report_unit >= 20 && (i) % report_unit == 0)
			show_msg(__func__, "Progress 2/2: %d/%d...\n", i, paths->len);
		path_i = g_ptr_array_index(paths, i); // If current path has similar seq with another alive path, remove current one.
		if (path_i->alive) {
			for (j = 0; j < paths->len; j++) {
				path_j = g_ptr_array_index(paths, j);
				if (path_i != path_j && path_j->alive) {
					similarity_score = similar_seqs(path_i->seq, path_j->seq,
							PATH_MISMATCHES * 4, PATH_MAX_GAPS, SCORE_MATCH,
							SCORE_MISMATCH, SCORE_GAP);
					if (similarity_score > 0) {
						// If similar, keep the longer one
						if (path_i->len < path_j->len)
							path_i->alive = 0;
						else
							path_j->alive = 0;
						break;
					}
				}
			}
		}
	}
	for (i = 0; i < paths->len; i++) {
		path_i = g_ptr_array_index(paths, i);
		if (!path_i->alive) {
			g_ptr_array_remove_fast(paths, path_i);
			i--;
		} else {
			//p_path(path_i);
		}
	}
}

GPtrArray *report_paths(edgearray *all_edges, bwa_seq_t *seqs) {
	int i = 0, j = 0;
	GPtrArray *block = NULL, *all_paths = NULL, *block_paths = NULL;
	edge *eg = NULL;

	all_paths = g_ptr_array_sized_new(1024);
	for (i = 0; i < all_edges->len; i++) {
		eg = g_ptr_array_index(all_edges, i);
		if (!eg->visited) {
			show_debug_msg(__func__,
					"REPORTING PATHS ---------------------------\n");
			block = g_ptr_array_sized_new(32);
			get_block_edges(eg, block);
			mark_duplicate_edges(block);
			show_debug_msg(__func__, "NEW BLOCK ---------------------------\n");
			for (j = 0; j < block->len; j++) {
				eg = g_ptr_array_index(block, j);
				p_flat_eg(eg);
			}
			if (get_level_n(block) >= MAX_ROADMAP_LEVEL) {
				show_debug_msg(__func__,
						"Too many level in this block (>=8) \n");
				block_paths = get_single_block_paths(block);
			} else {
				block_paths = iterate_block(block);
			}
			mark_duplicate_paths(block_paths);
			g_ptr_array_concat(all_paths, block_paths);
			g_ptr_array_free(block_paths, TRUE);
			g_ptr_array_free(block, TRUE);
		}
	}
	show_msg(__func__, "%d paths reported. \n", all_paths->len);
	determine_paths_ori(all_paths, seqs);
	mark_duplicate_paths(all_paths);
	show_msg(__func__, "%d paths after removing duplicates. \n", all_paths->len);
	return all_paths;
}

/**
 * Count how many commas in a string.
 */
int count_comma(const char *str, const int len) {
	int commas = 0, i = 0;
	for (i = 0; i < len; i++) {
		if (str[i] == ',')
			commas++;
	}
	return commas;
}

/**
 * Load the roadmap from txt files
 */
edgearray *load_rm(const hash_table *ht, const char *rm_dump_file,
		const char *rm_reads_file, const char *contig_file) {
	edgearray *edges = NULL;
	FILE *dump_fp = NULL, *reads_fp = NULL;
	int i = 0, j = 0, len = 0, no_right_connect_edge = 1;
	uint32_t n_ctgs = 0;
	int char_space = BUFSIZ - 1, char_len = 0;
	edge *eg = NULL, *in_out_eg = NULL, *eg_i = NULL;
	bwa_seq_t *r = NULL, *contigs = NULL, *ctg = NULL;
	char *read_str = (char*) calloc(char_space, sizeof(char)); // allocate buffer.
	char *attr[16], *shifts[3], ch = 0;

	dump_fp = xopen(rm_dump_file, "r");
	reads_fp = xopen(rm_reads_file, "r");
	edges = g_ptr_array_sized_new(BUFSIZ);
	show_msg(__func__, "Loading roadmap from %s and %s... \n", rm_dump_file,
			rm_reads_file);

	// Reconstruct the edges, pinpoint the reads on them first.
	while ((ch = fgetc(reads_fp)) != EOF) {
		if (char_len >= char_space) { // time to expand ?
			char_space += 1;
			kroundup32(char_space);
			// expand to double the current size of anything similar.
			read_str = realloc(read_str, char_space); // re allocate memory.
		}
		if (ch == '\n') {
			// Every line is in format of: [edge_id]	read_id_1,read_id_2,read_id_3...
			i = 0;
			attr[0] = strtok(read_str, "\t");
			while (attr[i] != NULL) { //ensure a pointer was found
				attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
			}
			eg = new_eg();
			eg->id = atoi(attr[0]);
			// show_debug_msg(__func__, "Initiating edge %s ...\n", attr[0]);

			if (attr[1] != NULL && strcmp(attr[1], "") != 0) {
				len = count_comma(attr[1], strlen(attr[1]));
				char *read_ids[len];
				read_ids[0] = strtok(attr[1], ",");
				i = 0;
				while (read_ids[i] != NULL) { //ensure a pointer was found
					read_ids[++i] = strtok(NULL, ","); //continue to tokenize the string
					//show_debug_msg(__func__, "%d: %s\n", i, read_ids[i]);
					if (read_ids[i] == NULL || strcmp(read_ids[i], "") == 0
							|| strcmp(read_ids[i], "\n") == 0) {
						break;
					}
					r = &ht->seqs[atoi(read_ids[i])];
					g_ptr_array_add(eg->reads, r);
					//p_query(__func__, r);
				}
			}
			g_ptr_array_add(edges, eg);
			free(read_str);
			read_str = (char*) calloc(BUFSIZ, sizeof(char)); // allocate buffer.
			char_space = BUFSIZ;
			char_len = 0;
		} else {
			read_str[char_len] = ch; // stuff in buffer.
			char_len++;
		}
	}
	show_msg(__func__, "Set up the connection among the edges... \n");
	free(read_str);
	read_str = (char*) calloc(BUFSIZ, sizeof(char)); // allocate buffer.
	char_space = BUFSIZ;
	char_len = 0;
	// Set up the connection among the edges
	for (i = 0; i < edges->len; i++) {
		eg = g_ptr_array_index(edges, i);
		while ((ch = fgetc(dump_fp)) != EOF) {
			if (char_len >= char_space) { // time to expand ?
				char_space += 1;
				kroundup32(char_space);
				// expand to double the current size of anything similar.
				read_str = realloc(read_str, char_space); // re allocate memory.
			}
			if (ch == '\n') {

				j = 0;
				attr[0] = strtok(read_str, "\t");
				while (attr[j] != NULL) { //ensure a pointer was found
					attr[++j] = strtok(NULL, "\t"); //continue to tokenize the string
				}

				// Outgoing edges
				if (attr[1] != NULL && strcmp(attr[1], "-1") != 0) {
					len = count_comma(attr[1], strlen(attr[1]));
					char *out_ids[len];
					out_ids[0] = strtok(attr[1], ",");
					j = 0;
					in_out_eg = edgearray_find_id(edges, atoi(out_ids[0]));
					if (in_out_eg)
						g_ptr_array_add(eg->out_egs, in_out_eg);
					while (out_ids[j] != NULL) { //ensure a pointer was found
						out_ids[++j] = strtok(NULL, ","); //continue to tokenize the string
						if (out_ids[j] == NULL || strcmp(out_ids[j], "") == 0) {
							break;
						}
						in_out_eg = edgearray_find_id(edges, atoi(out_ids[j]));
						if (in_out_eg)
							g_ptr_array_add(eg->out_egs, in_out_eg);
					}
				}
				// Incoming edges
				if (attr[2] != NULL && strcmp(attr[2], "-1") != 0) {
					len = count_comma(attr[2], strlen(attr[2]));
					char *in_ids[len];
					in_ids[0] = strtok(attr[2], ",");
					j = 0;
					in_out_eg = edgearray_find_id(edges, atoi(in_ids[0]));
					if (in_out_eg)
						g_ptr_array_add(eg->in_egs, in_out_eg);
					while (in_ids[j] != NULL) { //ensure a pointer was found
						in_ids[++j] = strtok(NULL, ","); //continue to tokenize the string
						if (in_ids[j] == NULL || strcmp(in_ids[j], "") == 0) {
							break;
						}
						in_out_eg = edgearray_find_id(edges, atoi(in_ids[j]));
						if (in_out_eg)
							g_ptr_array_add(eg->in_egs, in_out_eg);
					}
				}

				// Right contig
				if (attr[3] != NULL && strcmp(attr[3], "-1") != 0) {
					len = count_comma(attr[3], strlen(attr[3]));
					shifts[0] = strtok(attr[3], ","); //continue to tokenize the string
					shifts[1] = strtok(NULL, ",");
					show_debug_msg(__func__, "%s: %s\n", shifts[0], shifts[1]);
					in_out_eg = edgearray_find_id(edges, atoi(shifts[0]));
					if (in_out_eg) {
						eg->right_ctg = in_out_eg;
						g_ptr_array_free(eg->out_egs, TRUE);
						eg->out_egs = in_out_eg->out_egs;
						eg->r_shift = atoi(shifts[1]);
					}
				}
				free(read_str);
				read_str = (char*) calloc(BUFSIZ, sizeof(char)); // allocate buffer.
				char_space = BUFSIZ;
				char_len = 0;

				break; // Go to outside for loop to try next edge
			} else {
				read_str[char_len] = ch; // stuff in buffer.
				char_len++;
			}
		}
	}
	show_msg(__func__, "Loading the edge sequences... \n");
	// Assign the contig sequences
	contigs = load_reads(contig_file, &n_ctgs);
	for (i = 0; i < edges->len; i++) {
		no_right_connect_edge = 1;
		eg = g_ptr_array_index(edges, i);
		for (j = 0; j < n_ctgs; j++) {
			ctg = &contigs[j];
			if (eg->id == atoi(ctg->name)) {
				eg->contig = ctg;
				eg->len = ctg->len;
				break;
			}
		}
		if (eg->len <= 100) {
			g_ptr_array_remove_index_fast(edges, i);
			i--;
			continue;
		}
		// Set the root edges
		// If an edge right connects to some other edge, go to the edge
		while (eg->right_ctg) {
			eg = eg->right_ctg;
		}
		if (eg && eg->alive) {
			// If an edge has no incoming edges, set the be root
			if (eg->in_egs->len == 0)
				eg->is_root = 1;
			else {
				// Check all its input edges, if all of them right connect to this edge, set it to root
				for (j = 0; j < eg->in_egs->len; j++) {
					eg_i = g_ptr_array_index(eg->in_egs, j);
					if (!eg_i->right_ctg || eg_i->right_ctg != eg) {
						no_right_connect_edge = 0;
						break;
					}
				}
				if (no_right_connect_edge)
					eg->is_root = 1;
			}
		}
		//		p_ctg_seq("CONTIG", eg->contig);
	}
//	show_msg(__func__, "Updating reads... \n");
//	for (i = 0; i < edges->len; i++) {
//		eg = g_ptr_array_index(edges, i);
//		if (i % 100 == 0) {
//			show_msg(__func__, "Progress: %d/%d...\n", i, edges->len);
//		}
//		upd_reads_by_ol(ht->seqs, eg, MISMATCHES);
//	}
	free(read_str);
	fclose(reads_fp);
	fclose(dump_fp);
	return edges;
}

void align_back(hash_table *ht, rm_path *path, const int rl) {
	bwa_seq_t *query = NULL, *seqs = NULL, *read = NULL, *contig = NULL;
	alignarray *aligns = NULL;
	alg *a = NULL;
	int i = 0, j = 0, index = 0;

	seqs = ht->seqs;
	query = new_seq(path->seq, rl, 0);
	contig = path->seq;
	//p_query(__func__, query);
	aligns = g_ptr_array_sized_new(512);
	for (i = 0; i < path->len - rl; i++) {
		pe_aln_query(query, query->seq, ht, 2, query->len, 0, aligns);
		pe_aln_query(query, query->rseq, ht, 2, query->len, 1, aligns);
		for (j = 0; j < aligns->len; j++) {
			a = g_ptr_array_index(aligns, j);
			if (a->pos == 0) {
				index = a->r_id;
				read = &seqs[index];
				read->shift = i;
				read->rev_com = a->rev_comp;
				//p_query("ADDED", read);
				g_ptr_array_add(path->reads, read);
			}
		}
		reset_alg(aligns);
		ext_que(query, contig->seq[i], 0);
	}
	free_alg(aligns);
	bwa_free_read_seq(1, query);
}

/**
 * Stripe the starting and ending 0's of the path.
 *
 * If all coverage value is 0, return 0, indicating this path is not valid (no paired reads covering it)
 */
int trim_path(rm_path *path, int *coverage) {
	int i = 0, start_index = 0, end_index = path->len;
	bwa_seq_t *contig = path->seq;
	for (i = 0; i < path->len; i++) {
		if (coverage[i] == 0 && start_index == i) {
			start_index++;
		}
		if (coverage[i] != 0) {
			end_index = i + 1;
		}
	}
	if (start_index >= end_index)
		return 0;
	path->len = end_index - start_index;
	trun_seq(path->seq, start_index);
	contig->len = path->len;
	contig->seq[contig->len] = '\0';
	return 1;
}

/**
 * Assumption: the starting base and ending base with coverage larger than 0
 * For those cases with coverage pattern: ....000....000....
 * Break the path into several subpaths.
 */
void break_path(GPtrArray *paths, rm_path *path, const int *coverage) {
	int i = 0, flag = 1; // This flag indicates whether it is right time to break.
	rm_path *sub_path = NULL;
	int *pos_to_break = (int*) calloc(path->len, sizeof(int)), pos_no = 1;
	// Get the positions to break the path, in the format of (1_start, 1_end, 2_start, 2_end,...)
	for (i = 0; i < path->len; i++) {
		if (flag) {
			if (coverage[i] == 0) {
				pos_to_break[pos_no++] = i;
				flag = 0;
			}
		} else {
			if (coverage[i] > 0)
				flag = 1;
		}
	}
	if (pos_no % 2 != 0)
		pos_to_break[pos_no++] = path->len;
	if (pos_no > 2) {
		show_debug_msg(__func__, "Breaking this path: \n");
		//p_path(path);
		g_ptr_array_remove_fast(paths, path);
		for (i = 0; i < pos_no; i++) {
			if (i % 2 == 0) {
				sub_path = new_path();
				sub_path->id = path->id;
				// pos_to_break[i + 1] - pos_to_break[i] is the length
				sub_path->seq = new_seq(path->seq, pos_to_break[i + 1]
						- pos_to_break[i], pos_to_break[i]);
				sub_path->len = sub_path->seq->len;
				g_ptr_array_add(paths, sub_path);
			}
		}
	}
	free(pos_to_break);
}

int validate_p(hash_table *ht, GPtrArray *paths, const rm_path *path) {
	int i = 0, j = 0;
	char item[BUFSIZ];
	int n_forward = 0, n_backward = 0, n_counter_pairs = 0, n_reads = 0;
	int range_lower = 0, range_upper = 0;
	bwa_seq_t *left = NULL, *mate = NULL;
	int *coverage = NULL;
	readarray *left_mates = NULL;
	FILE *tmp = NULL;

	//	sprintf(item, "log/coverage_%d_%d.csv", path->id, path->len);
	//	tmp = xopen(item, "w");

	n_reads = path->reads->len;
	left_mates = g_ptr_array_sized_new(n_reads);
	coverage = (int*) calloc(path->len + 1, sizeof(int));
	for (i = 0; i < n_reads; i++) {
		left = g_ptr_array_index(path->reads, i);
		if (is_left_mate(left->name)) {
			g_ptr_array_add(left_mates, left);
			left->status = path->id;
		}
	}
	for (i = 0; i < n_reads; i++) {
		mate = g_ptr_array_index(path->reads, i);
		if (is_right_mate(mate->name)) {
			left = get_mate(mate, ht->seqs);
			if (left->status == path->id) {
				if (left->rev_com != mate->rev_com) {
					n_counter_pairs++;
				} else {
					if ((left->rev_com && left->shift >= mate->shift)
							|| (!left->rev_com && left->shift <= mate->shift)) {
						// The smaller value to be lower range
						if (left->shift <= mate->shift) {
							range_lower = left->shift;
							range_upper = mate->shift + mate->len;
						} else {
							range_lower = mate->shift;
							range_upper = left->shift + left->len;
						}
						for (j = range_lower; j < range_upper; j++) {
							coverage[j] += 1;
						}
						n_forward++;
					} else
						n_backward++;
				}
			}
		}
	}
	//	for (i = 0; i < path->len; i++) {
	//		sprintf(item, "%d\t%d\n", i, coverage[i]);
	//		fputs(item, tmp);
	//	}
	for (i = 0; i < left_mates->len; i++) {
		left = g_ptr_array_index(left_mates, i);
		left->status = 0;
	}
	trim_path(path, coverage); // Stripe the path
	break_path(paths, path, coverage); // Break the path if some point in the middle is 0
	//	fclose(tmp);
	g_ptr_array_free(left_mates, TRUE);
	free(coverage);
	return 1;
}

void validate_paths(hash_table *ht, GPtrArray *paths, const int rl) {
	int i = 0, is_valid = 1, report_unit = 0;
	rm_path *path;

	report_unit = paths->len / 1000;
	show_msg(__func__, "Validating paths... \n");
	for (i = 0; i < paths->len; i++) {
		if (report_unit != 0 && (i) % report_unit == 0)
			show_msg(__func__, "Progress: %d/%d...\n", i, paths->len);
		path = g_ptr_array_index(paths, i);
		align_back(ht, path, rl);
		is_valid = validate_p(ht, paths, path);
		if (!is_valid) {
			path->alive = 0;
			g_ptr_array_remove_index_fast(paths, i);
			i--;
		}
	}
}

int test_sw(const char *fa_fn) {
	bwa_seq_t *seqs = NULL;
	int n_seqs = 0, score = 0;
	seqs = load_reads(fa_fn, &n_seqs);
	p_ctg_seq("SEQ1", &seqs[0]);
	p_ctg_seq("SEQ1", &seqs[1]);
	score = similar_seqs(&seqs[0], &seqs[1], 16, 4, 2, -1, -1);
	show_debug_msg(__func__, "Score: %d \n", score);
	return 1;
}

int pe_path(int argc, char *argv[]) {
	clock_t t = clock();
	hash_table *ht = NULL;
	edgearray *edges = NULL;
	GPtrArray *final_paths = NULL;
	FILE *merged_pair_contigs = NULL;

	show_msg(__func__, "%s \n", argv[1]);
	show_msg(__func__, "%s \n", argv[2]);
	show_msg(__func__, "%s \n", argv[3]);
	show_msg(__func__, "%s \n", argv[4]);
	show_msg(__func__, "%s \n", argv[5]);

	if (!g_thread_supported())
		g_thread_init(NULL);
	ht = pe_load_hash(argv[2]);
	edges = load_rm(ht, argv[3], argv[4], argv[5]);
	reset_edge_ids(edges);
	merge_ol_edges(edges, 197, 50, ht, 1);
	merged_pair_contigs
			= xopen("../SRR097897_out/pair_contigs.1.fa", "w");
	save_edges(edges, merged_pair_contigs, 0, 0, 100);
	fclose(merged_pair_contigs);
	//scaffolding(edges, 197, ht, 4);
	//graph_by_edges(edges, "../SRR097897_part/roadmap.1.dot");
	//dump_rm(edges, "../SRR097897/roadmap.graph", "../SRR097897/roadmap.reads");
	//final_paths = report_paths(edges, ht->seqs);
	//validate_paths(ht, final_paths, atoi(argv[1]));
	//save_paths(final_paths, "../SRR097897_out/peta.fa", 100);
	//g_ptr_array_free(final_paths, TRUE);

	// test_sw(argv[4]);

	fprintf(stderr, "[pe_path] Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
