/*
 * hash.h
 *
 *  Created on: 09-May-2013
 *      Author: carl
 */

#ifndef HASH_H_
#define HASH_H_

#include <unordered_map>
#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwtaln.h"
#include "utils.h"
#include "k_hash.h"
#include "tpl.hpp"

using namespace std;
typedef unordered_map<uint64_t, uint64_t*> kmer_hash;
typedef unordered_map<int, tpl*> tpl_hash;
typedef unordered_map<uint64_t, uint64_t> mer_counter;

typedef struct {
	uint64_t kmer;
	uint64_t count;
} kmer_counter;

#ifdef __cplusplus
extern "C" {
#endif

	void build_tpl_hash(kmer_hash &hash, tpl_hash *tpls, const int k,
			const int read_len);
	gint cmp_kmers_by_count(gpointer a, gpointer b);

#ifdef __cplusplus
}
#endif

#endif /* HASH_H_ */
