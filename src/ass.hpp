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
#include "tpl.hpp"
#include "hash.hpp"

using namespace std;

typedef struct {
	hash_table *ht;
	tpl_hash *all_tpls;
	int to_try_connect;
} kmer_t_meta;


#ifdef __cplusplus
extern "C" {
#endif

	int pe_kmer(int argc, char *argv[]);
	int pe_cluster(int argc, char *argv[]);
	void strip_branches(hash_table *ht, tpl_hash *all_tpls, tpl *t);
	void finalize_tpl(hash_table *ht, tpl_hash *all_tpls, tpl *t, int to_branching, int to_con_left,
			int to_con_right);

#ifdef __cplusplus
}
#endif

#endif /* ASS_HPP_ */
