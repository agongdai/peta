/*
 * pool.h
 *
 *  Created on: 20-Jun-2011
 *      Author: carl
 */
#ifndef POOL_H_
#define POOL_H_
#include <glib.h>
#include "bwtaln.h"
#include "edge.h"

#ifdef __cplusplus
extern "C" {
#endif
#define POOLSIZE 32

typedef struct {
	readarray *reads;
	size_t n;
	size_t space;
} pool;

pool *new_pool();
void free_pool(pool *r_pool);
void p_edge_reads(const edge *eg);
void p_readarray(const readarray *ra, const int part);
void p_pool(const char *header, const pool *r_pool, const int *next);
int binary_exists(const pool *r_pool, const bwa_seq_t *read);
int get_next_char(pool *cur_pool, const int ori, edge *eg);
int pool_exists(const pool *p, const bwa_seq_t *read);
void pool_add(pool *p, bwa_seq_t *new_seq);
void mate_pool_add(pool *p, bwa_seq_t *new_seq);
void pool_uni_add(pool *p, bwa_seq_t *new_seq);
bwa_seq_t *forward(pool *cur_pool, const char c, edge *ass_eg, const int left_max_ctg_id);
void pool_sort_ins(pool *r_pool, bwa_seq_t *new_seq);
void rm_partial(pool *cur_pool, int ori, bwa_seq_t *query, int nm);
void clean_mate_pool(pool *mate_pool);
gboolean pool_rm(pool *r_pool, bwa_seq_t *rm_seq);
gboolean pool_rm_fast(pool *p, bwa_seq_t *read);
gboolean pool_rm_index(pool *p, const int i);
void syn_pools(pool *cur_pool, pool *mate_pool, const bwa_seq_t *seqs, const int ori);
void clear_pool(pool *r_pool);
void pool_get_majority(pool *cur_pool, const char c, edge *ass_eg);

#ifdef __cplusplus
}
#endif

#endif /* POOL_H_ */
