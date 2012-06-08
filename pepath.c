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
#include "pepath.h"
#include "roadmap.h"
#include "bwase.h"
#include "pechar.h"

int path_id = 0;

extern unsigned char nst_nt4_table[256];

void p_path(const rm_path *p) {
	int i = 0;
	edge *contig;
	node *n;
	bwa_seq_t *query;
	GList *eg_item = 0;
	if (!p) {
		printf("Empty path! \n");
		return;
	}
	printf("[p_path] ----------------------------------------\n");
	printf("[p_path] Path %d (%p): %d \n", p->id, p, p->len);
	eg_item = g_list_first(p->edges);
	while (eg_item && i++ <= p->n_ctgs) {
		contig = eg_item->data;
		if (!contig) {
			printf("[p_path] Contig null in path %d \n", p->id);
			continue;
		}
		n = contig->left;
		query = contig->query;
		if (i == 1 || i == p->n_ctgs)
			printf("[p_path] \t Contig %d [%s]: %d (%d, %d)\n", contig->id,
					query->name, contig->len,
					contig->right_ctg ? contig->right_ctg->id : 0,
					contig->r_shift);
		else
			printf("[p_path] \t Contig %d: %d (%d, %d)\n", contig->id,
					contig->len, contig->right_ctg ? contig->right_ctg->id : 0,
					contig->r_shift);
		eg_item = eg_item->next;
	}
	printf("[p_path] ----------------------------------------\n");
}

void save_path(const rm_path *path, FILE *ass_fa) {
	bwa_seq_t *contig;
	edge *eg;
	GList *eg_item = 0;
	int i = 0, j = 0, start_pos = 0, shift = 0;
	char *h = malloc(BUFSIZE);
	char char_n = to_upper_lower('N'), c = char_n;
	if (!path || !ass_fa)
		return;
	eg_item = g_list_first(path->edges);
	eg = eg_item->data;
	printf("[save_path] Writing path %d [%s] to disk... \n", path->id,
			eg->query->name);
	sprintf(h, ">%d|len=%d \n", path->id, path->len);
	//	fprintf(stderr, "[save_path] Writing path %d (%p) to disk... %d\n", path->id, path,
	//			path->len);
	fputs(h, ass_fa);
	while (eg_item && i++ <= path->n_ctgs) {
		eg = eg_item->data;
		contig = eg->contig;
		if (!contig || !contig->seq) {
			fprintf(
					stderr,
					"[save_path] Contig sequence null here for edge %d, query %s. \n",
					eg->id, eg->query->name);
			break;
		}

		for (j = shift; j < eg->len; j++) {
			if (j < 0) {
				fputc(char_n, ass_fa);
			} else {
				c = contig->seq[j];
				if (nst_nt4_table[(int) c] <= 4 && nst_nt4_table[(int) c] >= 0) {
					// c = contig->seq[j];
				} else
					c = char_n;
				fputc(c, ass_fa);
			}
			if ((start_pos + j + 1 - shift) % LINELEN == 0) {
				fputc('\n', ass_fa);
			}
		}
		start_pos = (start_pos + j + 1 - shift) % LINELEN - 1;
		shift = eg->right_ctg ? eg->r_shift : 0;
		eg_item = eg_item->next;
	}
	if ((i - 1) % LINELEN != 0)
		fputc('\n', ass_fa);
	free(h);
}

void upd_path_len(rm_path *path) {
	int i = 0, total_len = 0, shift = 0;
	edge *eg, *next_eg;
	GList *eg_item = 0, *eg_item_next = 0;
	if (!path)
		return;
	eg_item = g_list_first(path->edges);
	while (eg_item && i++ < path->n_ctgs) {
		eg = eg_item->data;
		//		printf("[%d] eg->right_ctg = %p, eg->shift = %d. \n", path->id,
		//				eg->right_ctg, eg->shift);
		// If current edge is the last edge from left roadmap,
		//	its next edge should be the first edge from right roadmap,
		//	whose right_ctg is current one, so here set the length value.
		if (i < path->n_ctgs - 1) {
			eg_item_next = eg_item->next;
			if (eg_item_next) {
				next_eg = eg_item_next->data;
				if (next_eg && next_eg->right_ctg && next_eg->right_ctg->id
						== eg->id) {
					eg->len = next_eg->r_shift;
					next_eg->right_ctg = 0;
					next_eg->r_shift = 0;
				}
			}
		}
		total_len += eg->len - shift;
		if (eg->right_ctg) {
			shift = eg->r_shift;
		} else {
			shift = 0;
		}
		eg_item = eg_item->next;
		//		printf("[%d] eg->len = %d, shift = %d, total_len = %d. \n", path->id,
		//				eg->len, shift, total_len);
	}
	path->len = total_len;
}

void save_paths(const rm_path *paths, const int n_paths, FILE *ass_fa) {
	rm_path *p;
	int i = 0;
	if (!paths || !ass_fa)
		return;
	fprintf(stderr, "[save_paths] Writing %d paths to disk... \n", n_paths);
	for (i = 0; i < n_paths; i++) {
		p = (rm_path*) &paths[i];
		if (p) {
			upd_path_len(p);
			p_path(p);
			save_path(p, ass_fa);
		}
	}
}

rm_path *new_path() {
	rm_path *p = (rm_path*) malloc(sizeof(rm_path));
	p->edges = (edgelist*) calloc(1, sizeof(edgelist));
	p->len = 0;
	p->n_ctgs = 0;
	p->id = ++path_id;
	return p;
}

rm_path *dup_path(rm_path *p) {
	rm_path *p2 = (rm_path*) malloc(sizeof(rm_path));
	if (!p)
		return 0;
	p2->edges = g_list_copy(p->edges);
	p2->len = p->len;
	p2->id = ++path_id;
	p2->n_ctgs = p->n_ctgs;
	return p2;
}

void rvs_path(rm_path *p) {
	int n = 0, i = 0;
	edge *eg;
	bwa_seq_t *contig;
	GList *eg_item = 0;
	if (!p)
		return;
	n = p->n_ctgs;
	p->edges = g_list_reverse(p->edges);
	eg_item = g_list_first(p->edges);
	while (eg_item && i++ <= n) {
		eg = eg_item->data;
		contig = eg->contig;
		seq_reverse(contig->len, contig->seq, 0);
		eg_item = eg_item->next;
	}
}

int app_ctg_to_path(rm_path *p, edge *contig) {
	if (!contig)
		return 1;
	if (!p)
		p = new_path();
	if (g_list_find(p->edges, contig)) {
		return 1;
	}

	p->edges = g_list_append(p->edges, contig);
	p->len += contig->len;
	p->n_ctgs++;
	return 0;
}

void add_path(rm_path *all, rm_path *p, int *n_paths, int *space) {
	int n = *n_paths;
	int p_index = 0;
	rm_path *new_paths, *tmp;
	if (!all || !p)
		return;
	// If the path id exists, ignore it.
	for (p_index = 0; p_index < n; p_index++) {
		tmp = &all[p_index];
		if (tmp->id == p->id)
			return;
	}
	// If there is no enough space, double the space.
	if (path_id >= *space - 1) {
		new_paths = (rm_path*) malloc(sizeof(rm_path) * (*space * 2));
		memcpy(new_paths, all, sizeof(rm_path) * (*space));
		*all = *new_paths;
		*space *= 2;
	}
	all[n] = *p;
	*n_paths += 1;
}

void trc_path(rm_path *p, const int pos) {
	int new_len = 0;
	edge *eg;
	GList *eg_item = 0;
	if (!p || pos >= p->n_ctgs || pos < 0)
		return;
	eg_item = g_list_first(p->edges);
	while (eg_item) {
		eg = eg_item->data;
		if (eg)
			new_len += eg->len;
		eg_item = eg_item->next;
	}
	p->n_ctgs = pos + 1;
	eg_item = g_list_nth(p->edges, pos + 1);
	eg_item->next = 0;
	p->len = new_len;
}

void merge_paths(rm_path *left_path, const int pos, const rm_path *right_path) {
	if (!left_path || !right_path || pos >= left_path->n_ctgs)
		return;
	trc_path(left_path, pos);
	left_path->edges = g_list_concat(left_path->edges, right_path->edges);
	left_path->len += right_path->len;
	left_path->n_ctgs += right_path->n_ctgs;
}

void rtv_edge_paths(edge *eg, rm_path *path, rm_path *all_paths, int *n_paths,
		int *space, int left_max_ctg_id) {
	node *right;
	rm_path *p2;
	edge *tmp, *right_eg;
	int i = 0, repeated = 0;
	if (!eg || !all_paths)
		return;
	right = eg->right;
	repeated = app_ctg_to_path(path, eg);
	if (repeated)
		return;
	right_eg = eg->right_ctg;
	if (right_eg) {
		repeated = app_ctg_to_path(path, right_eg);
		if (repeated)
			return;
		path->len -= eg->r_shift;
		if (eg->id > left_max_ctg_id && right_eg->id <= left_max_ctg_id)
			return;
	}
	for (i = 0; i < right->out; i++) {
		tmp = g_ptr_array_index(right->edge_out, i);
		if (tmp && tmp->id <= left_max_ctg_id && eg->id > left_max_ctg_id)
			return;
		if (i == right->out - 1) {
			rtv_edge_paths(tmp, path, all_paths, n_paths, space,
					left_max_ctg_id);
		} else {
			p2 = dup_path(path);
			rtv_edge_paths(tmp, p2, all_paths, n_paths, space, left_max_ctg_id);
			add_path(all_paths, p2, n_paths, space);
		}
	}
}

rm_path *rtv_rm_paths(const roadmap *rm, int *n_paths,
		const int left_max_ctg_id) {
	rm_path *paths, *p, *p2;
	node *node;
	edge *start_e, *eg, *right_eg;
	int i = 0, j = 0, space = INIT_PATH_N, repeated = 0;
	if (!rm)
		return 0;
	paths = (rm_path*) calloc(space, sizeof(rm_path));
	start_e = rm->start_egs;
	for (i = 0; i < rm->start_eg_n; i++) {
		start_e = g_ptr_array_index(rm->start_egs, i);
		p = new_path();
		app_ctg_to_path(p, start_e);
		// The 'right_ctg' is just partially on the path, cut by property 'cursor'.
		right_eg = start_e->right_ctg;

		if (right_eg) {
			if (start_e->id > left_max_ctg_id && right_eg->id
					<= left_max_ctg_id) {
				add_path(paths, p, n_paths, &space);
				continue;
			}
			repeated = app_ctg_to_path(p, right_eg);
			if (repeated)
				continue;
			p->len -= start_e->r_shift;
		}
		node = start_e->right;
		for (j = 0; j < node->out; j++) {
			eg = g_ptr_array_index(node->edge_out, j);
			if (j == node->out - 1) {
				rtv_edge_paths(eg, p, paths, n_paths, &space, left_max_ctg_id);
			} else {
				p2 = dup_path(p);
				rtv_edge_paths(eg, p2, paths, n_paths, &space, left_max_ctg_id);
				// Adding path must be done AFTER the path is completed!
				add_path(paths, p2, n_paths, &space);
			}
		}
		// Not put in the for loop, in case there is no out-going edges at all
		add_path(paths, p, n_paths, &space);
	}
	return paths;
}

int is_pet_wrapped(rm_path *path) {
	edge *first_eg, *last_eg;
	char *left_query_name, *right_query_name;
	bwa_seq_t *first_seq, *last_seq;
	if (!path)
		return 0;
	first_eg = g_list_first(path->edges)->data;
	last_eg = g_list_last(path->edges)->data;
	if (!first_eg || !last_eg)
		return 0;
	first_seq = first_eg->contig;
	last_seq = last_eg->contig;
	if (!first_seq || !last_seq)
		return 0;
	left_query_name = first_seq->name;
	right_query_name = last_seq->name;
	//	printf("[is_pet_wrapped] Path %d: [%s, %s] \n", path->id, left_query_name,
	//			right_query_name);
	return is_mates(left_query_name, right_query_name);
}

void rtv_back_paths(edge *eg, rm_path *path, rm_path *all_paths, int *n_paths,
		int *space, int level) {
	node *left;
	rm_path *p2;
	edge *tmp;
	int i = 0, repeated = 0;
	if (!eg || !all_paths)
		return;
	if (level-- <= 0)
		return;
	left = eg->right;
	repeated = app_ctg_to_path(path, eg);
	if (repeated)
		return;
	for (i = 0; i < left->in; i++) {
		tmp = g_ptr_array_index(left->edge_in, i);
		if (i == left->in - 1) {
			rtv_back_paths(tmp, path, all_paths, n_paths, space, level);
		} else {
			p2 = dup_path(path);
			rtv_back_paths(tmp, p2, all_paths, n_paths, space, level);
			add_path(all_paths, p2, n_paths, space);
		}
	}
}

rm_path *mer_lr_paths(rm_path *left_paths, rm_path *right_paths,
		const int left_n_paths, const int right_n_paths, int *n_paths) {
	int i = 0, j = 0, k = 0, space = INIT_PATH_N, shift = 0;
	rm_path *paths, *p, *p_right, *p2;
	edge *eg, *r_eg, *last_eg, *sc_last_eg;
	GList *eg_item = 0;
	if (!left_paths || !right_paths)
		return 0;

	// Reverse all right paths, as well as their contigs;
	for (i = 0; i < right_n_paths; i++) {
		p = &right_paths[i];
		if (p) {
			rvs_path(p);
		}
	}

	paths = (rm_path*) calloc(space, sizeof(rm_path));
	// @TODO: Three nested for loop, could be improved.
	for (i = 0; i < left_n_paths; i++) {
		p = &left_paths[i];
		// p_path(p);
		if (p) {
			// An edge from right roadmap could connect to any edge in the left roadmap
			eg_item = g_list_first(p->edges);
			while (eg_item && j++ <= p->n_ctgs) {
				eg = eg_item->data;
				// Check if there is any path from right roadmap connects here
				for (k = 0; k < right_n_paths; k++) {
					p_right = &right_paths[k];
					// Since it has been reversed, the first edge contains the right_ctg info
					r_eg = p_right->edges;
					if (r_eg && r_eg->right_ctg && r_eg->right_ctg->id
							== eg->id) {
						// A left path could be used multiple times
						p2 = dup_path(p);
						last_eg = &p2->edges[p2->n_ctgs - 1];
						// r_eg->right_ctg = 0; // To avoid infinite loop, one right path could be used once.
						if (last_eg) {
							// Update the path length
							p2->len -= (last_eg->len - r_eg->r_shift);
							last_eg->contig->len = shift;
						}
						// Check the second last edge;
						// If its right contig is the last edge, compare the shift value against
						//	the right roadmap edge. If its shift value is large, abandon it.
						if (p2->n_ctgs > 1) {
							sc_last_eg = &p2->edges[p2->n_ctgs - 2];
							if (sc_last_eg && sc_last_eg->right_ctg
									&& sc_last_eg->r_shift > r_eg->r_shift)
								continue;
						}
						merge_paths(p2, j - 1, p_right);
						if (is_pet_wrapped(p2))
							add_path(paths, p2, n_paths, &space);
					}
				}
				eg_item = eg_item->next;
			}
		}
	}
	return paths;
}

rm_path *rtv_paths(const roadmap *left_rm, const roadmap *right_rm,
		const int left_max_ctg_id, int *n_paths) {
	rm_path *left_paths, *right_paths, *paths, *p;
	int left_n_paths = 0, right_n_paths = 0, i = 0;

	if (!left_rm || !left_rm || left_max_ctg_id <= 0)
		return 0;
	left_paths = rtv_rm_paths(left_rm, &left_n_paths, left_max_ctg_id);
	right_paths = rtv_rm_paths(right_rm, &right_n_paths, left_max_ctg_id);
	printf("[rtv_paths] %d left paths \n", left_n_paths);
	printf("[rtv_paths] %d right paths \n", right_n_paths);

	for (i = 0; i < left_n_paths; i++) {
		p = &left_paths[i];
		p_path(p);
	}

	paths = mer_lr_paths(left_paths, right_paths, left_n_paths, right_n_paths,
			n_paths);
	return paths;
}

