/*
 * graph.hpp
 *
 *  Created on: Jun 30, 2013
 *      Author: carl
 */

#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <glib.h>
#include "tpl.h"
#include "bwtaln.h"
#include "hash.hpp"

typedef struct {
	int id;					// Vertex id
	bwa_seq_t *ctg;			// Sequences.
	int len;				// Length of vertex, same as ctg->len
	float weight;			// Likely to be # of reads
	GPtrArray *ins;			// Incoming edges
	GPtrArray *outs;		// Outgoing edges
	GPtrArray *reads;		// Reads on it, for paired-end tracking
	int status;				// 0 means good.
	int index;				// To determine strongly connected components in the graph
	int lowlink;			// http://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
} vertex;

typedef struct {
	int id;					// Edge id
	bwa_seq_t *junc_seq;	// Sequence in the junction: max read_len * 2 - 8
	int left_len;			// The junc_seq is from left vertex and right vertex,
	int right_len;			//		whose lengths are not always the same
	int len;				// Length of the edge
	vertex *left;			// Left vertex
	vertex *right;			// Right vertex
	float weight;			// Likely to be # of reads
	GPtrArray *reads;		// Reads on it, for paired-end tracking
	int status;				// 0 means good.
} edge;

typedef struct {
	GPtrArray *vertexes;
	GPtrArray *edges;
	int len;
	int index;				// Global index for Tarjan's algorithm
	GPtrArray *stack;		// Stack for Tarjan's algorithm
	GPtrArray *scc;			// Strongly connected components
} splice_graph;

typedef struct {
	GPtrArray *vertexes;
	GPtrArray *edges;
} component;

typedef struct {
	vertex *pre;
	vertex *post;
	GPtrArray *vertexes;
	GPtrArray *edges;
} ASM;

void assign_reads2tpl(tpl *eg, hash_map *hm);
void process_graph(GPtrArray *all_tpls, GPtrArray *all_juncs,
		hash_map *hm);
void p_vertex(vertex *v);
GPtrArray *reads_on_seq(bwa_seq_t *seq, hash_map *hm, const int n_mismatch);
int vertex_in_scc(splice_graph *g, vertex *v);

#endif /* GRAPH_HPP_ */
