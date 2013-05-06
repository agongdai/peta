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
#define		MIN_BRANCH_LEN		8

using namespace std;
typedef unordered_map<uint64_t, edge*> tpl_hash;

typedef struct {
	hash_map *hm;
	tpl_hash *all_tpls;
} kmer_t_meta;

typedef struct {
	edge *main_tpl;
	edge *branch_tpl;
	int locus;
	int weight;
	uint8_t ori;
	uint64_t kmer;	// Kmer when branching
} junction;

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
	void kmer_ext_edge(edge *eg, uint64_t query_int, hash_map *hm, tpl_hash *all_tpls, const int ori);

#ifdef __cplusplus
}
#endif

#endif /* ASS_HPP_ */
