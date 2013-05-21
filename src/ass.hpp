/*
 * ass.hpp
 *
 *  Created on: 09-Apr-2013
 *      Author: carl
 */

#ifndef ASS_HPP_
#define ASS_HPP_

#include <glib.h>
#include <unordered_map>
#include "kmers.hpp"
#include "edge.h"

#define		SHORT_BRANCH_SHIFT	4
#define 	MIN_WEIGHT			2
#define		ERROR_RATE			0.02
#define		BRANCH_THRE			0.1

using namespace std;
typedef unordered_map<int, edge*> tpl_hash;

typedef struct {
	hash_map *hm;
	tpl_hash *all_tpls;
} kmer_t_meta;

typedef struct {
	edge *tpl;
	int locus;
	int8_t main_c;
	int main_count;
	int8_t second_c;
	int second_count;
} ambi_base;

#ifdef __cplusplus
extern "C" {
#endif

	int pe_kmer(int argc, char *argv[]);
	int kmer_ext_edge(edge *eg, uint64_t query_int, hash_map *hm, tpl_hash *all_tpls, const int ori);

#ifdef __cplusplus
}
#endif

#endif /* ASS_HPP_ */
