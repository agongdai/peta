#ifndef PATH_HPP_
#define PATH_HPP_

#include <glib.h>
#include <unordered_map>
#include "tpl.hpp"
#include "bwtaln.h"
#include "hash.hpp"
#include "graph.hpp"
#include "k_hash.h"

using namespace std;
typedef unordered_map<vertex*, float> vertex_p_hash;

typedef struct {
	int id;					// Path id
	GPtrArray *vertexes;	// Vertexes
	GPtrArray *edges;		// Edges
	GPtrArray *reads;		// Reads on it
	bwa_seq_t *ctg;			// Contig sequence
	int len;				// Length of path
	int *junction_points;	// Junction points where reads validation should be done
	int *junction_lengths;	// Junction lengths
	float *weights;			// Weights of vertexes and junctions, same junction may have different weights on paths
	float coverage;			// Coverage
	uint8_t status;			// 0 means good
} path;

void determine_paths(splice_graph *g, hash_table *ht, char *save_dir);

#endif /* GRAPH_HPP_ */
