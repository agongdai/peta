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
} vertex;

typedef struct {
	int id;					// Edge id
	bwa_seq_t *junc_seq;	// Sequence in the junction: max read_len * 2 - 8
	int left_len;			// The junc_seq is from left vertex and right vertex,
	int right_len;			//		whose lengths are not always the same
	vertex *left;			// Left vertex
	vertex *right;			// Right vertex
	float weight;			// Likely to be # of reads
	GPtrArray *reads;		// Reads on it, for paired-end tracking
} edge;

typedef struct {
	GPtrArray *vertexes;	// Vertexes
	GPtrArray *edges;		// Edges
	bwa_seq_t *ctg;			// Contig sequence
	int len;				// Length of path
} path;

typedef struct {
	GPtrArray *vertexes;
	GPtrArray *edges;
} splice_graph;

void assign_reads2tpl(tpl *eg, hash_map *hm);
void process_graph(GPtrArray *all_tpls, GPtrArray *all_juncs,
		hash_map *hm);

#endif /* GRAPH_HPP_ */
