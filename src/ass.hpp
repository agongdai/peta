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
#include "tpl.hpp"
#include "read.h"

using namespace std;
typedef unordered_map<int, tpl*> tpl_hash;

typedef struct {
	hash_map *hm;
	tpl_hash *all_tpls;
} kmer_t_meta;

#ifdef __cplusplus
extern "C" {
#endif

	int pe_kmer(int argc, char *argv[]);
	int kmer_ext_tpl(tpl *t, uint64_t query_int, hash_map *hm, tpl_hash *all_tpls, const int ori);

#ifdef __cplusplus
}
#endif

#endif /* ASS_HPP_ */
