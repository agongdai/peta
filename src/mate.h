/*
 * mate.h
 *
 *  Created on: 18-Mar-2013
 *      Author: carl
 */

#ifndef MATE_H_
#define MATE_H_

#include "pehash.h"
#include "pool.h"
#include "edgelist.h"
#include "bwtaln.h"

#define N_BIG_MATE_POOL			10000000
#define SD_TIMES 			 	4
#define MAX_COUNTER_PAIRS		4

#ifdef __cplusplus
extern "C" {
#endif

	void keep_pairs_only(edge *eg, bwa_seq_t *seqs);
	void add_mates_by_ol(const bwa_seq_t *seqs, edge *eg, pool *cur_pool,
			const int ol, const int nm, bwa_seq_t *query, const int ori,
			const int insert_size, const int sd_insert_size);
	pool *get_mate_pool_from_edge(edge *eg, const bwa_seq_t *seqs, const int ori,
			const int insert_size, const int sd_insert_size);

#ifdef __cplusplus
}
#endif

#endif /* MATE_H_ */
