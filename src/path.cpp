#include <stdint.h>
#include <glib.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "rnaseq.h"
#include "pechar.h"
#include "utils.h"
#include "peseq.h"
#include "bwtaln.h"
#include "tpl.hpp"
#include "kmers.hpp"
#include "hash.hpp"
#include "graph.hpp"
#include "path.hpp"
#include "junction.hpp"
#include "k_hash.h"

int path_id = 0;

path *new_path() {
	path *p = (path*) malloc(sizeof(path));
	p->id = path_id++;
	p->vertexes = g_ptr_array_sized_new(0);
	p->edges = g_ptr_array_sized_new(0);
	p->reads = NULL;
	p->len = 0;
	p->ctg = NULL;
	p->junction_points = NULL;
	p->weights = NULL;
	p->junction_lengths = NULL;
	p->coverage = 0;
	p->status = 0;
	return p;
}

path *clone_path(path *p) {
	path *clone = new_path();
	uint32_t i = 0;
	vertex *v = NULL;
	edge *e = NULL;
	float *weights = NULL, w = 0.0;
	for (i = 0; i < p->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(p->vertexes, i);
		g_ptr_array_add(clone->vertexes, v);
	}
	for (i = 0; i < p->edges->len; i++) {
		e = (edge*) g_ptr_array_index(p->edges, i);
		g_ptr_array_add(clone->edges, e);
	}
	return clone;
}

void destroy_path(path *p) {
	if (p) {
		if (p->edges)
			g_ptr_array_free(p->edges, TRUE);
		if (p->vertexes)
			g_ptr_array_free(p->vertexes, TRUE);
		if (p->reads)
			g_ptr_array_free(p->reads, TRUE);
		if (p->ctg)
			bwa_free_read_seq(1, p->ctg);
		if (p->junction_points)
			free(p->junction_points);
		if (p->weights)
			free(p->weights);
		if (p->junction_lengths)
			free(p->junction_lengths);
		free(p);
	}
}

void destroy_paths(GPtrArray *paths) {
	uint32_t i = 0;
	path *p = NULL;
	if (paths) {
		for (i = 0; i < paths->len; i++) {
			p = (path*) g_ptr_array_index(paths, i);
			destroy_path(p);
		}
		g_ptr_array_free(paths, TRUE);
	}
}

void destory_levels(GPtrArray *levels) {
	GPtrArray *vertexes = NULL;
	uint32_t i = 0;
	if (levels) {
		for (i = 0; i < levels->len; i++) {
			vertexes = (GPtrArray*) g_ptr_array_index(levels, i);
			g_ptr_array_free(vertexes, TRUE);
		}
		g_ptr_array_free(levels, TRUE);
	}
}

void p_p(path *p) {
	uint32_t i = 0;
	vertex *v = NULL;
	edge *e = NULL;
	show_debug_msg(__func__, "==== Path %d ====\n", p->id);
	printf("\tLength: %d\n", p->len);
	printf("\tVertexes: %d \n", p->vertexes->len);
	for (i = 0; i < p->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(p->vertexes, i);
		printf("\t\tVertex %d: %d; Weight %.2f\n", v->id, v->len, v->weight);
	}
	printf("\tEdges: %d \n", p->edges->len);
	for (i = 0; i < p->edges->len; i++) {
		e = (edge*) g_ptr_array_index(p->edges, i);
		printf("\t\tEdge %d: %d; Weight %.2f\n", e->id, e->len, e->weight);
	}
	if (p->weights) {
		printf("\tWeights: %d \n", p->edges->len + p->vertexes->len);
		for (i = 0; i < p->edges->len + p->vertexes->len; i++) {
			printf("\t\tWeight: %.2f \n", p->weights[i]);
		}
	}
}

void p_paths(GPtrArray *paths, char *fn) {
	path *p = NULL;
	uint32_t i = 0, j = 0;
	vertex *v = NULL;
	edge *e = NULL;
	char entry[BUFSIZ];
	FILE *dot = xopen(fn, "w");
	fputs("digraph G {\n", dot);
	fputs("graph [rankdir=LR];\n", dot);
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		if (p->status == 0)
		sprintf(
				entry,
				"%d [label=\"Path %d\" style=filled fillcolor=\"yellow\" shape=box]; \n",
				p->id, p->id);
		else
			sprintf(
					entry,
					"%d [label=\"Path %d\" style=filled fillcolor=\"gray\" shape=box]; \n",
					p->id, p->id);
		fputs(entry, dot);
		for (j = 0; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			if (j == 0) {
				sprintf(entry, "%d -> %d.%d; \n", p->id, v->id, i);
				fputs(entry, dot);
			}
			sprintf(entry, "%d.%d [label=\"%d: %d\" shape=box]; \n", v->id, i,
					v->id, v->len);
			fputs(entry, dot);
		}
		for (j = 0; j < p->edges->len; j++) {
			e = (edge*) g_ptr_array_index(p->edges, j);
			sprintf(entry, "%d.%d -> %d.%d [label=\"%.0f\"]; \n", e->left->id,
					i, e->right->id, i, p->weights[2 * j + 1]);
			fputs(entry, dot);
		}
	}
	fputs("}\n", dot);
	fclose(dot);
}

void save_paths(GPtrArray *paths, char *fn, const int to_save_all) {
	uint32_t i = 0;
	path *p = NULL;
	char entry[BUFSIZ];
	FILE *p_fp = xopen(fn, "w");
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		// If we are not saving all paths and the status is not 0, skip
		if (!to_save_all) {
            if (p->status != 0 || p->len < 100)
			    continue;
        }
		sprintf(entry, ">%d length: %d\n", p->id, p->len);
		save_con(entry, p->ctg, p_fp);
	}
	fclose(p_fp);
}

/*
 * Get junction sequences
 */
bwa_seq_t *get_breaking_seq(path *p, int breaking_index, int read_len) {
	bwa_seq_t *seq = p->ctg;
	int *points = p->junction_points;
	int n_points = p->edges->len;
	int half_max = read_len - SHORT_BRANCH_SHIFT;
	int start = points[breaking_index] - half_max, end = points[breaking_index]
			+ half_max;
	int shared_start = 0, shared_end = 0, shared_size;
	int tuned_end = 0, tuned_start = 0;
	// Make sure not to double count any read between adjacent junctions
	//show_debug_msg(__func__, "Start: %d; End: %d\n", start, end);
	if (n_points > breaking_index + 1) {
		// Size of shared region of current junction and next junction (if any)
		shared_size = end - (points[breaking_index + 1] - half_max);
		if (shared_size >= read_len) {
			shared_start = points[breaking_index + 1] - half_max;
			shared_end = shared_start + shared_size;
			tuned_end = ((shared_end - (read_len - 1)) + shared_start) / 2
					+ (read_len - 1);
			end = tuned_end < end ? tuned_end : end;
			/**
			 show_debug_msg(__func__, "Shared: [%d -> %d]\n", shared_start,
			 shared_end);
			 show_debug_msg(__func__, "Tuned end at breaking index %d: %d \n",
			 breaking_index, end);
			 **/
		}
	}
	if (breaking_index > 0) {
		shared_size = points[breaking_index - 1] + half_max - start;
		if (shared_size >= read_len) {
			shared_start = start;
			shared_end = shared_start + shared_size;
			tuned_start = ((shared_end - (read_len - 1)) + shared_start) / 2;
			start = tuned_start > start ? tuned_start : start;
			/**
			 show_debug_msg(__func__, "Shared: [%d -> %d]\n", shared_start,
			 shared_end);
			 show_debug_msg(__func__, "Tuned start at breaking index %d: %d \n",
			 breaking_index, start);
			 **/
		}
	}
	start = (start < 0) ? 0 : start;
	end = (end >= seq->len) ? seq->len : end;
	//show_debug_msg(__func__, "Start: %d; End: %d\n", start, end);
	return new_seq(seq, end - start, start);
}

/**
 * Get vertex without incoming edges
 */
GPtrArray *get_graph_roots(comp *c) {
	vertex *v = NULL;
	uint32_t i = 0;
	GPtrArray *roots = g_ptr_array_sized_new(16);
	for (i = 0; i < c->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(c->vertexes, i);
		if (!(v->ins) || v->ins->len == 0)
			g_ptr_array_add(roots, v);
	}
	return roots;
}

/**
 * Get levels of vertexes.
 * Return an array of array:
 * 	[0]: vertexes at level 0
 *  [1]: vertexes at level 1
 *  ......
 */
GPtrArray *get_vertex_levels(comp *c) {
	GPtrArray *level_vertexes = get_graph_roots(c), *next_level = NULL;
	GPtrArray *levels = g_ptr_array_sized_new(32);
	uint32_t i = 0, j = 0, has_more = 1, n_level = 0;
	vertex *v = NULL;
	edge *e = NULL;
	/**
	 for (i = 0; i < g->vertexes->len; i++) {
	 v = (vertex*) g_ptr_array_index(g->vertexes, i);
	 show_debug_msg(__func__, "Vertex [%d, %d] in graph \n", v->id, v->len);
	 }
	 **/
	for (i = 0; i < level_vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(level_vertexes, i);
		v->status = 0;
	}
	while (has_more) {
		g_ptr_array_add(levels, level_vertexes);
		has_more = 0;
		next_level = NULL;
		n_level++;
		show_debug_msg(__func__, "Level %d: %d \n", n_level, level_vertexes->len);
		for (i = 0; i < level_vertexes->len; i++) {
			v = (vertex*) g_ptr_array_index(level_vertexes, i);
			//p_vertex(v);
			if (v->outs->len > 0) {
				if (!next_level)
					next_level = g_ptr_array_sized_new(level_vertexes->len);
				for (j = 0; j < v->outs->len; j++) {
					e = (edge*) g_ptr_array_index(v->outs, j);
					//p_edge(e);
					// If the edge is not visited before
					// If the vertex is added before at the same level, ignore
					if (e->right->status != n_level) {
						g_ptr_array_add(next_level, e->right);
						e->right->status = n_level;
						has_more = 1;
					}
				}
			}
		}
		level_vertexes = next_level;
	}
	return levels;
}

/**
 * Get combinatorial paths from levels of vertexes.
 */
GPtrArray *combinatorial_paths(GPtrArray *levels) {
	GPtrArray *level_vertexes = NULL;
	path *p = NULL, *new_p = NULL;
	GPtrArray *paths = g_ptr_array_sized_new(levels->len), *next_paths = NULL;
	uint32_t i = 0, j = 0, k = 0, m = 0;
	vertex *v = NULL, *tail_v = NULL;
	edge *e = NULL;
	if (levels->len == 0)
		return paths;
	level_vertexes = (GPtrArray*) g_ptr_array_index(levels, 0);
	// Create
	for (j = 0; j < level_vertexes->len; j++) {
		v = (vertex*) g_ptr_array_index(level_vertexes, j);
		p = new_path();
		g_ptr_array_add(p->vertexes, v);
		g_ptr_array_add(paths, p);
	}
	// Append the vertexes to current paths level by level
	for (i = 1; i < levels->len; i++) {
		level_vertexes = (GPtrArray*) g_ptr_array_index(levels, i);
		next_paths = g_ptr_array_sized_new(paths->len * 2);
		// For paths with last vertex connected to this vertex
		for (k = 0; k < paths->len; k++) {
			p = (path*) g_ptr_array_index(paths, k);
			new_p = NULL;
			tail_v = (vertex*) g_ptr_array_index(p->vertexes, p->vertexes->len
					- 1);
			// For each vertex at current level
			for (j = 0; j < level_vertexes->len; j++) {
				v = (vertex*) g_ptr_array_index(level_vertexes, j);
				// Check whether the last vertex of the path is connected to the current vertex
				for (m = 0; m < tail_v->outs->len; m++) {
					e = (edge*) g_ptr_array_index(tail_v->outs, m);
					if (e->right == v) {
						new_p = clone_path(p);
						g_ptr_array_add(new_p->vertexes, v);
						g_ptr_array_add(new_p->edges, e);
						g_ptr_array_add(next_paths, new_p);
					}
				}
			}
			// Since path p is freed later, make a clone and add
			if (!new_p)
				g_ptr_array_add(next_paths, clone_path(p));
		}
		// Free the paths from last level, assign next level paths to current one
		destroy_paths(paths);
		paths = next_paths;

		/**
		 show_debug_msg(__func__, "Paths after level %d\n", i);
		 for (j = 0; j < paths->len; j++) {
		 p = (path*) g_ptr_array_index(paths, j);
		 p_p(p);
		 }
		 **/
	}
	show_debug_msg(__func__, "%d Paths reported \n", paths->len);
	return paths;
}

/**
 * If an exon is shorter than k,
 * 		verify whether there are reads covering this vertex and its left and right vertexes
 */
void validate_short_exons(GPtrArray *paths, hash_table *ht) {
	uint32_t i = 0, j = 0;
	int start = 0, end = 0, n_reads = 0;
	bwa_seq_t *seq = NULL;
	GPtrArray *reads = NULL;
	int read_len = ht->o->read_len;
	path *p = NULL;
	vertex *v = NULL;
	edge *e = NULL;
	float w = 0.0;
	// Destroy path with any edge weight 0
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		//p_p(p);
		for (j = 0; j < p->edges->len; j++) {
			w = p->weights[j * 2 + 1];
			if (w <= 0) {
				show_debug_msg(__func__, "Path %d is not valid with 0 edge \n",
						p->id);
				destroy_path(p);
				g_ptr_array_remove_index_fast(paths, i--);
				break;
			}
		}
	}
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		start = end = 0;
		if (p->vertexes->len < 2)
			continue;
		// Ignore the first exon
		v = (vertex*) g_ptr_array_index(p->vertexes, 0);
		start = v->len;
		for (j = 1; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			// If the vertex is smaller than the k value (25)
			if (v->len < ht->o->k + N_MISMATCHES) {
				end = start + read_len - N_MISMATCHES * 2;
				end = end > p->len ? p->len : end;
				start = start - (read_len - v->len) + N_MISMATCHES * 2;
				start = start < 0 ? 0 : start;
				seq = new_seq(p->ctg, end - start, start);
				reads = reads_on_seq(seq, ht, N_MISMATCHES);

				/**
				 show_debug_msg(__func__, "Path %d exon %d: [%d, %d)\n", p->id, v->id, start, end);
				 p_ctg_seq("PATH", p->ctg);
				 p_ctg_seq("EXON", seq);
				 p_readarray(reads, 0);
				**/

				n_reads = reads->len;
				g_ptr_array_free(reads, TRUE);
				bwa_free_read_seq(1, seq);
				if (n_reads == 0) {
					show_debug_msg(__func__, "Path %d is not valid \n", p->id);
					destroy_path(p);
					g_ptr_array_remove_index_fast(paths, i--);
					break;
				}
			}
			start += v->len;
		}
	}
}

/**
 * Initialize attributes to paths
 */
void assign_path_attrs(GPtrArray *paths, hash_table *ht) {
	uint32_t i = 0, j = 0, len = 0;
	path *p = NULL;
	vertex *v = NULL;
	int *points = NULL;
	edge *e = NULL;
	bwa_seq_t *seq = NULL;
	GPtrArray *reads = NULL;
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		//show_debug_msg(__func__, "Properties of path %d \n", p->id);
		p->len = 0;
		// Get path length
		for (j = 0; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			p->len += v->len;
		}
		// Concatenate the vertexes to path
		p->ctg = blank_seq(p->len);
		len = 0;
		for (j = 0; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			memcpy(p->ctg->seq + len, v->ctg->seq, sizeof(ubyte_t) * v->len);
			len += v->len;
		}
		p->ctg->len = len;
		p->len = len;
		// The coverage, given all reads are from this path
		p->reads = reads_on_seq(p->ctg, ht, N_MISMATCHES);
		g_ptr_array_sort(p->reads, (GCompareFunc) cmp_reads_by_name);
		p->coverage = ((float) p->reads->len) * ((float) ht->o->read_len)
				/ ((float) p->len);
		// Determine the junction points
		points = (int*) malloc(sizeof(int) * (p->vertexes->len - 1));
		len = 0;
		for (j = 0; j < p->vertexes->len - 1; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			len += v->len;
			points[j] = len;
		}
		p->junction_points = points;
		// Assign the weights of vertexes and edges
		p->weights = (float*) calloc((p->edges->len + p->vertexes->len),
				sizeof(float));
		p->junction_lengths = (int*) calloc(p->edges->len, sizeof(int));
		for (j = 0; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			p->weights[2 * j] = (float) v->weight;
			//show_debug_msg(__func__, "Weight of vertex %d: %.2f\n", v->id,
			//		v->weight);
		}
		for (j = 0; j < p->edges->len; j++) {
			e = (edge*) g_ptr_array_index(p->edges, j);
			seq = get_breaking_seq(p, j, ht->o->read_len);
			// Stringent in junctions, 0 mismatches.
			reads = reads_on_seq(seq, ht, N_MISMATCHES);

			/**
			 show_debug_msg(__func__, "=== Path %d, Breaking point: %d ===\n",
			 p->id, points[j]);
			 p_ctg_seq("PATH", p->ctg);
			 p_ctg_seq(__func__, seq);
			 p_readarray(reads, 1);
			 show_debug_msg(__func__, "=== %d reads. END === \n\n", reads->len);
			 **/

			p->weights[2 * j + 1] = (float) reads->len;
			//show_debug_msg(__func__, "Weight of edge %d: %.2f\n", e->id,
			//		p->weights[2 * j + 1]);
			p->junction_lengths[j] = seq->len;
			e->len = seq->len;
			g_ptr_array_free(reads, TRUE);
			bwa_free_read_seq(1, seq);
		}
	}
}

/**
 * Return initial path probability;
 */
float *init_path_prob(GPtrArray *paths, const float read_len) {
	float cov_sum = 0.0;
	path *p = NULL;
	uint32_t i = 0, j = 0;
	vertex *v = NULL;
	float *p_cov = (float*) calloc(paths->len, sizeof(float));
	float *init_path_p = (float*) calloc(paths->len, sizeof(float));
	if (paths->len == 0) {
		free(p_cov);
		return init_path_p;
	}
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		p_cov[i] = 0.0;
		for (j = 0; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			p_cov[i] += (v->weight * (float) read_len) / ((float) v->len);
		}
		for (j = 0; j < p->edges->len; j++) {
			p_cov[i] += p->weights[2 * j + 1] * read_len
					/ p->junction_lengths[j];
		}
		cov_sum += p_cov[i];
	}
	if (cov_sum == 0) {
		for (i = 0; i < paths->len; i++) {
			init_path_p[i] = 1.0 / (float) paths->len;
		}
	} else {
		for (i = 0; i < paths->len; i++) {
			p = (path*) g_ptr_array_index(paths, i);
			init_path_p[i] = p_cov[i] / cov_sum;
		}
	}
	free(p_cov);
	return init_path_p;
}

int array_contains(GPtrArray *arr, gpointer element) {
	uint32_t i = 0;
	gpointer p = NULL;
	for (i = 0; i < arr->len; i++) {
		p = g_ptr_array_index(arr, i);
		if (p == element)
			return 1;
	}
	return 0;
}

/**
 * Input: 'cvoerage' an array of float arrays:
 * Array:
 * 		Path 1: cov(feature_1_0), cov(feature_1_1)...
 * 		Path 2: cov(feature_2_0), cov(feature_2_1)...
 * 		...
 * Every path contains: exon_1, junction_1, exon_2, junction_2, ... exon_n, junction_n, exon_n+1
 * 		n_exons = n_junctions + 1
 */
GPtrArray *calc_features_coverage(GPtrArray *paths, const int read_len,
		float *path_p, GPtrArray *coverage) {
	float *cov = NULL, f_cov = 0.0, f_p_sum = 0.0, weight = 0.0;
	uint32_t i = 0, j = 0, m = 0;
	path *p = NULL, *p2 = NULL;
	vertex *v = NULL;
	edge *e = NULL;
	if (!coverage) {
		// If the value is NULL, means initializing the coverage array
		coverage = g_ptr_array_sized_new(paths->len);
		for (i = 0; i < paths->len; i++) {
			p = (path*) g_ptr_array_index(paths, i);
			cov = (float*) calloc(p->vertexes->len + p->edges->len,
					sizeof(float));
			g_ptr_array_add(coverage, cov);
		}
	}
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		cov = (float*) g_ptr_array_index(coverage, i);
		// Calc the coverage of every vertex on this path
		for (j = 0; j < p->vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(p->vertexes, j);
			weight = v->weight > 0 ? v->weight : 0;
			f_cov = weight * read_len / v->len;
			f_p_sum = 0;
			// Sum path prob of all paths containing this vertex
			for (m = 0; m < paths->len; m++) {
				p2 = (path*) g_ptr_array_index(paths, m);
				if (array_contains(p2->vertexes, (gpointer) v)) {
					f_p_sum += path_p[m];
					//show_debug_msg(__func__, "Path %d contains exon %d; f_p_sum = %.2f\n", p2->id, v->id, f_p_sum);
				}
			}
			//show_debug_msg(__func__, "Path %d, exon %d: %.2f \n", p->id, v->id, f_p_sum);
			cov[j * 2] = f_cov * path_p[i] / f_p_sum;
		}
		// Calc the coverage of every edge on this path
		// The weights are stored in p->weights
		f_p_sum = 0;
		//p_p(p);
		for (j = 0; j < p->edges->len; j++) {
			//show_debug_msg(__func__, "Path %d, edge %d \n", p->id, j);
			e = (edge*) g_ptr_array_index(p->edges, j);
			weight = p->weights[j * 2 + 1];
			weight = weight > 0 ? weight : 0;
			f_cov = weight * read_len / p->junction_lengths[j];
			f_p_sum = 0;
			for (m = 0; m < paths->len; m++) {
				p2 = (path*) g_ptr_array_index(paths, m);
				if (array_contains(p2->edges, (gpointer) e))
					f_p_sum += path_p[m];
			}
			cov[j * 2 + 1] = f_cov * path_p[i] / f_p_sum;
		}
	}
	return coverage;
}

/**
 * Parameters:
 * 		paths_p: Probability of paths
 * 		g: splice graph
 * 		paths: combinatorial paths
 * 		read_len: read length
 */
void diffsplice_em(comp *c, GPtrArray *paths, const float read_len,
		float *paths_p) {
	int round = 0, n_features = 0, max_iteration = 200000, stop = 0;
	float stop_threshold = 0.000001;
	uint32_t i = 0, j = 0, m = 0, n = 0, n_paths = paths->len;
	// Coverage of paths
	float *p_covs = NULL;
	// Feature coverage on a transcript
	float *cov = NULL;
	GPtrArray *coverage = NULL;
	// Total # of reads on the the paths
	float sum_N = 0.0, sum_p_len = 0.0, *this_Ns = (float*) calloc(
			sizeof(float), n_paths), sum_p = 0.0, sum_len = 0;
	float *next_paths_p = NULL, sum_next_f_p = 0.0, sum_next_f_t_p = 0.0;
	float *k_te = NULL, sum_k_te = 0.0, sum_k_c_te = 0.0, sum_denominator = 0.0;

	path *p = NULL, *p2 = NULL;
	vertex *v = NULL;

	if (paths->len <= 1)
		return;

	p_covs = (float*) calloc(sizeof(float), n_paths);
	// Hash: vertex => vertex_len/sum_len
	vertex_p_hash v_p_h;
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		sum_p_len += p->len;
		for (j = 0; j < p->edges->len; j++) {
			sum_p_len += p->junction_lengths[j];
		}
	}

	//show_debug_msg(__func__, "Total length: %.2f\n", sum_p_len);
	// Set the vertex probability: vertex_len/sum_len
	for (i = 0; i < c->vertexes->len; i++) {
		v = (vertex*) g_ptr_array_index(c->vertexes, i);
		v_p_h[v] = ((float) v->len) / ((float) sum_p_len);
		//printf("Vertex %d: %d/%.2f, %.2f \n", v->id, v->len, sum_p_len, v_p_h[v]);
	}

	// Set initial coverage values
	coverage = calc_features_coverage(paths, read_len, paths_p, coverage);

	show_debug_msg(__func__, "========== Initialization ==========\n");
	for (i = 0; i < n_paths; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		show_debug_msg(__func__, "---- Path %d ---- \n", p->id);
		printf("\t\t\tProbability: %.2f \n", paths_p[i]);
	}

	while (round++ < max_iteration && !stop) {
		// Expectation step: E-step
		for (i = 0; i < n_paths; i++) {
			p = (path*) g_ptr_array_index(paths, i);
			//show_debug_msg(__func__, "---- Path %d ---- \n", p->id);
			n_features = p->vertexes->len + p->edges->len;
			k_te = (float*) calloc(n_features, sizeof(float));
			// Set k(t,e)
			for (j = 0; j < p->vertexes->len; j++) {
				v = (vertex*) g_ptr_array_index(p->vertexes, j);
				k_te[2 * j] = (p->len * v->len) / read_len * (p->len - v->len);
			}
			for (j = 0; j < p->edges->len; j++) {
				k_te[2 * j + 1] = (p->len * p->junction_lengths[j]) / read_len
						* (p->len - p->junction_lengths[j]);
			}

			for (j = 0; j < n_features; j++) {
				//show_debug_msg(__func__, "k_te[%d]: %.2f\n", j, k_te[j]);
			}

			sum_k_te = 0.0;
			for (j = 0; j < n_features; j++) {
				sum_k_te += k_te[j];
			}

			sum_k_c_te = 0.0;
			cov = (float*) g_ptr_array_index(coverage, i);
			for (j = 0; j < n_features; j++) {
				sum_k_c_te += k_te[j] * (cov[j] * cov[j]);
				//show_debug_msg(__func__, "cov[%d, %d]: %.2f\n", i, j, cov[j]);
			}

			//show_debug_msg(__func__, "sum_k_c_te[i]: %.2f\n", i, sum_k_c_te);

			// Coverage of path i in this iteration
			p_covs[i] = (0 - n_features + sqrt(n_features * n_features + 4
					* sum_k_te * sum_k_c_te)) / (2 * sum_k_te);

			//show_debug_msg(__func__, "p_covs[i]: %.2f\n", i, p_covs[i]);

			// # of reads from path i in this iteration
			sum_len = 0;
			for (j = 0; j < p->edges->len; j++) {
				sum_len += p->junction_lengths[j];
			}
			this_Ns[i] = p_covs[i] * (sum_len + p->len) / read_len;

			free(k_te);
		}

		/**
		 printf("\n\n==== After E-step %d ====\n", round);
		 for (i = 0; i < n_paths; i++) {
		 p = (path*) g_ptr_array_index(paths, i);
		 show_debug_msg(__func__, "---- Path %d ---- \n", p->id);
		 printf("\t\t\tCoverage: %.2f \n", p_covs[i]);
		 printf("\t\t\tN_reads: %.2f \n", this_Ns[i]);
		 printf("\t\t\tProbability: %.2f \n", paths_p[i]);
		 }
		 **/

		// Sum of reads from all paths
		sum_N = 0;
		for (j = 0; j < n_paths; j++) {
			sum_N += this_Ns[j];
		}

		//printf("N: %.2f\n", sum_N);

		// Maximization step: M-step
		next_paths_p = (float*) calloc(sizeof(float), n_paths);
		for (j = 0; j < n_paths; j++) {
			p = (path*) g_ptr_array_index(paths, j);
			sum_next_f_p = 0.0;
			for (m = 0; m < n_paths; m++) {
				p2 = (path*) g_ptr_array_index(paths, m);
				if (m == j)
					continue;
				sum_next_f_t_p = 0.0;
				for (n = 0; n < p2->vertexes->len; n++) {
					v = (vertex*) g_ptr_array_index(p2->vertexes, n);
					sum_next_f_t_p += v_p_h[v];
					//printf("sum_next_f_t_p[%d, %d, %d]: %.2f \n", j, m, n, sum_next_f_t_p);
				}
				for (n = 0; n < p2->edges->len; n++) {
					sum_next_f_t_p += ((float) p2->junction_lengths[n])
							/ ((float) sum_p_len);
				}
				sum_next_f_p += paths_p[m] * sum_next_f_t_p;
			}

			//show_debug_msg(__func__, "sum_next_f_p[%d]: %.2f\n", j, sum_next_f_p);

			sum_next_f_t_p = 0.0;
			for (n = 0; n < p->vertexes->len; n++) {
				v = (vertex*) g_ptr_array_index(p->vertexes, n);
				sum_next_f_t_p += v_p_h[v];
			}
			for (n = 0; n < p->edges->len; n++) {
				sum_next_f_t_p += ((float) p->junction_lengths[n])
						/ ((float) sum_p_len);
				//show_debug_msg(__func__,
				//		"junction_lengths[%d] / sum_p_len: %d/%d += %.2f\n", n,
				//		p->junction_lengths[n], sum_p_len, sum_next_f_t_p);
			}

			/**
			show_debug_msg(__func__, "sum_next_f_p[%d]: %.2f\n", j,
					sum_next_f_p);
			show_debug_msg(__func__, "sum_next_f_t_p[%d]: %.2f\n", j,
					sum_next_f_t_p);
			show_debug_msg(__func__, "this_Ns[%d]/Sum: %.2f / %.2f \n", j,
					this_Ns[j], sum_N);
			**/

			sum_denominator = (sum_next_f_t_p * (sum_N - this_Ns[j]));
			if (sum_denominator > 0.0) {
				next_paths_p[j] = (this_Ns[j] * sum_next_f_p) / sum_denominator;
			} else {
				if (this_Ns[j] >= sum_N)
					next_paths_p[j] = 1.0;
				else
					next_paths_p[j] = 0.0;
			}
		}
		// Sum and normalize the next path prob
		sum_p = 0.0;
		for (j = 0; j < n_paths; j++) {
			sum_p += next_paths_p[j];
		}
		for (j = 0; j < n_paths; j++) {
			next_paths_p[j] = next_paths_p[j] / sum_p;
		}
		stop = 1;
		for (j = 0; j < n_paths; j++) {
			/**
			show_debug_msg(__func__, "%.5f vs %.5f \n", paths_p[j],
					next_paths_p[j]);
			show_debug_msg(__func__, "Difference: %.5f \n", get_abs(paths_p[j]
					- next_paths_p[j]));
			**/
			if ((float) (get_abs(paths_p[j] - next_paths_p[j]))
					> stop_threshold) {
				stop = 0;
			}
		}
		if (!stop) {
			for (j = 0; j < n_paths; j++) {
				paths_p[j] = next_paths_p[j];
			}
		}
		free(next_paths_p);

		coverage = calc_features_coverage(paths, read_len, paths_p, coverage);

		/**
		printf("\n\n==== Iteration %d ====\n", round);
		for (i = 0; i < n_paths; i++) {
			p = (path*) g_ptr_array_index(paths, i);
			show_debug_msg(__func__, "---- %d: Path %d ---- \n", i, p->id);
			printf("\t\t\tCoverage: %.2f \n", p_covs[i]);
			printf("\t\t\tN_reads: %.2f \n", this_Ns[i]);
			printf("\t\t\tProbability: %.3f \n", paths_p[i]);
		}
		**/

	}
	show_debug_msg(__func__, "Stop at Iteration %d \n", round);
	for (i = 0; i < n_paths; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		show_debug_msg(__func__, "---- %d: Path %d ---- \n", i, p->id);
		printf("\t\t\tCoverage: %.2f \n", p_covs[i]);
		printf("\t\t\tN_reads: %.2f \n", this_Ns[i]);
		printf("\t\t\tProbability: %.3f \n", paths_p[i]);
	}
	for (i = 0; i < coverage->len; i++) {
		cov = (float*) g_ptr_array_index(coverage, i);
		free(cov);
	}
	g_ptr_array_free(coverage, TRUE);
	free(p_covs);
	free(this_Ns);
}

gint cmp_paths_by_cov(gpointer a, gpointer b) {
    path *p_a = *((path**) a);
    path *p_b = *((path**) b);
    if (p_a->coverage >= p_b->coverage)
        return 1;
    else
        return -1;
}

void high_cov_paths(GPtrArray *paths, const int max_n_paths) {
    uint32_t i = 0;
    path *p = NULL;
    g_ptr_array_sort(paths, (GCompareFunc) cmp_paths_by_cov);
    for (i = 0; i < paths->len; i++) {
        if (i >= max_n_paths - 1)
            break;
        p = (path*) g_ptr_array_index(paths, i);
        show_debug_msg(__func__, "Path [%d, %d] with coverage %.2f obtained\n", p->id, p->len, p->coverage);
    }
    while(paths->len > max_n_paths) 
        g_ptr_array_remove_index_fast(paths, paths->len - 1);
}

/**
 * Determine the paths of the component
 */
GPtrArray *comp_paths(comp *c, hash_table *ht) {
	show_debug_msg(__func__, "Getting vertexes in component %d...\n", c->id);
	GPtrArray *levels = get_vertex_levels(c);
	uint32_t j = 0, i = 0;
	vertex *v = NULL;
	GPtrArray *vertexes = NULL;
	float *paths_prob = NULL;
	path *p = NULL;
	for (i = 0; i < levels->len; i++) {
		vertexes = (GPtrArray*) g_ptr_array_index(levels, i);
		show_debug_msg(__func__, ">>>> Level %d <<<<<\n", i);
		for (j = 0; j < vertexes->len; j++) {
			v = (vertex*) g_ptr_array_index(vertexes, j);
			p_vertex(v);
		}
		printf("\n");
	}
	GPtrArray *paths = combinatorial_paths(levels);
	destory_levels(levels);
	assign_path_attrs(paths, ht);
	validate_short_exons(paths, ht);

    if (paths->len > 2000) {
        high_cov_paths(paths, 3);
        return paths;
    }

	paths_prob = init_path_prob(paths, ht->o->read_len);
	diffsplice_em(c, paths, ht->o->read_len, paths_prob);
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		// If the probability of the paths is small, mark it as not alive.
		if (paths_prob[i] <= 0.02 && paths_prob[i] * paths->len < 0.5)
			p->status = 1;
	}
	free(paths_prob);
	return paths;
}

void append_paths(GPtrArray *all_paths, GPtrArray *paths) {
	uint32_t i = 0;
	path *p = NULL;
	for (i = 0; i < paths->len; i++) {
		p = (path*) g_ptr_array_index(paths, i);
		g_ptr_array_add(all_paths, p);
	}
}

GPtrArray *paths_from_tpl(comp *c) {
    tpl *from = NULL;
    uint32_t i = 0;
    vertex *v = NULL;
    path *p = NULL;
    GPtrArray *paths = g_ptr_array_sized_new(4);
    for (i = 0; i < c->vertexes->len; i++) {
        v = (vertex*) g_ptr_array_index(c->vertexes, i);
        from = v->from;
        if (from->alive) {
            show_debug_msg(__func__, "Component too complicated, graph template [%d, %d] only. \n", from->id, from->len);
            p = new_path();
            p->ctg = new_seq(from->ctg, from->len, 0);
            g_ptr_array_add(paths, p);
            // Mark it as 'dead', in case duplicated
            from->alive = 0;
        }
    }
    return paths;
}

void determine_paths(splice_graph *g, hash_table *ht, char *save_dir) {
	comp *c = NULL;
	uint32_t i = 0;
	path *p = NULL;
	char *fn = NULL;
	GPtrArray *paths = NULL, *all_paths = g_ptr_array_sized_new(
			g->vertexes->len / 10);
	for (i = 0; i < g->components->len; i++) {
		c = (comp*) g_ptr_array_index(g->components, i);
        p_comp(c, save_dir);
		if (c->vertexes->len >= MAX_VS_IN_COMP)
            paths = paths_from_tpl(c);
        else 
		    paths = comp_paths(c, ht);
		append_paths(all_paths, paths);
	}
	fn = get_output_file("paths.fa", save_dir);
	save_paths(all_paths, fn, 1);
	free(fn);
	fn = get_output_file("peta.fa", save_dir);
	save_paths(all_paths, fn, 0);
	free(fn);

	fn = get_output_file("paths.dot", save_dir);
	p_paths(all_paths, fn);
	free(fn);
	destroy_paths(all_paths);
}
