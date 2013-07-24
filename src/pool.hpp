/*
 * pool.h
 *
 *  Created on: Jul 23, 2013
 *      Author: carl
 */

#ifndef POOL_H_
#define POOL_H_
#include <glib.h>
#include "bwtaln.h"
#include "k_hash.h"
#include "tpl.hpp"

using namespace std;

// A simple encapsulation on list of reads
typedef struct {
	GPtrArray *reads;
} pool;

#ifdef __cplusplus
extern "C" {
#endif

	pool *new_pool();
	void destroy_pool(pool *p);
	void p_pool(const char *header, const pool *r_pool, const int *next);

	void add2pool(pool *p, bwa_seq_t *r);
	int get_next_char(pool *p, tpl *t, const int ori);
	void forward(pool *p, tpl *t, const int ori);
	void next_pool(hash_table *ht, pool *p, tpl *t, bwa_seq_t *tail,
			int mismatches, const int ori);
	void correct_tpl_base(pool *p, tpl *t, int t_len);
	void rm_half_clip_reads(pool *p, tpl *t, int tpl_c, int mismatches, int ori);

#ifdef __cplusplus
}
#endif

#endif /* POOL_H_ */
