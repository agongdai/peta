/*
 * ass.hpp
 *
 *  Created on: 09-Apr-2013
 *      Author: carl
 */

#ifndef ASS_HPP_
#define ASS_HPP_

#include <glib.h>
#include "kmers.hpp"

typedef struct {
	hash_map *hm;
	GPtrArray *all_edges;
} kmer_t_meta;

int pe_kmer(int argc, char *argv[]);

#endif /* ASS_HPP_ */
