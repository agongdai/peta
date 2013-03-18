/*
 * mate.c
 *
 *  Created on: 18-Mar-2013
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "pool.h"
#include "pehash.h"
#include "bwase.h"
#include "utils.h"
#include "peseq.h"
#include "mate.h"
#include "edgelist.h"

/**
 * From current edge, get all mates of the used reads.
 */
pool *get_mate_pool_from_edge(edge *eg, const hash_table *ht, const int ori,
		const int insert_size, const int sd_insert_size) {
	int i = 0;
	bwa_seq_t *s = NULL, *mate = NULL, *seqs = NULL;
	pool *mate_pool = NULL;
	mate_pool = new_pool();
	seqs = ht->seqs;
	for (i = 0; i < eg->reads->len; i++) {
		s = g_ptr_array_index(eg->reads, i);
		mate = get_mate(s, seqs);
		//p_query("READ", s);
		//p_query("MATE", mate);
		// If the mate should have been used.
		if (is_paired(s, ori))
			continue;
		// If the insert size is not in the range
		if (abs(eg->len - s->shift) > (insert_size + sd_insert_size * SD_TIMES)) {
			continue;
		}
		// If the mate is already in use, either by current or another thread
		if (mate->is_in_c_pool || mate->is_in_m_pool || mate->status == USED
				|| mate->status == DEAD)
			continue;
		// If the used read is used by another thread;
		//	or the mate has been used by this template before.
		if (!(s->status == TRIED && s->contig_id == eg->id) || (mate->status
				== TRIED && mate->contig_id == eg->id))
			continue;
		mate->rev_com = s->rev_com;
		mate_pool_add(mate_pool, mate, eg->tid);
	}
	return mate_pool;
}

/**
 * Add read to current pool from the mate pool.
 * Check whether a mate overlaps with the tail with length parameter 'nm'
 */
void add_mates_by_ol(const hash_table *ht, edge *eg, pool *cur_pool,
		const int ol, const int nm, bwa_seq_t *query, const int ori,
		const int insert_size, const int sd_insert_size) {
	int i = 0, overlapped = 0;
	bwa_seq_t *mate = NULL, *tmp = NULL;
	readarray *ol_mates = NULL;
	bwa_seq_t *template = NULL, *seqs = NULL;
	pool *mate_pool = NULL;
	reads_ht *rht = NULL;
	seqs = ht->seqs;
	// Copy read length of the end of the contig.
	template = new_seq(eg->contig, seqs->len, eg->len - seqs->len);
	if (ori) {
		seq_reverse(template->len, template->seq, 0);
	}
	mate_pool = get_mate_pool_from_edge(eg, ht, ori, insert_size,
			sd_insert_size);
	//p_readarray(mate_pool->reads, 1);
	if (mate_pool->n >= N_BIG_MATE_POOL) {
		rht = build_reads_ht(ol, mate_pool->reads);
		ol_mates = find_reads_ol_template(rht, template, seqs, ori);
	} else {
		ol_mates = mate_pool->reads;
	}
	//p_readarray(ol_mates, 1);
	// Add the mate reads which overlap with the tail into the current pool
	for (i = 0; i < ol_mates->len; i++) {
		mate = g_ptr_array_index(ol_mates, i);
		// For the reads in the mate pool, these ones are not considered:
		//	1. Is already in c_pool (in this or another thread), or in the mate pool of other thread
		//	2. Is already used
		//	3. Its mate is not used
		//	4. Its mate is used, by by another edge
		//	5. The distance between the mates are out of range
		tmp = mate;
		if (mate->rev_com)
			tmp = new_mem_rev_seq(mate, mate->len, 0);
		overlapped = find_ol_within_k(tmp, template, nm, ol - 1,
				query->len - 1, ori);
		/*if (strcmp(mate->name, "2460877") == 0) {
		 show_debug_msg("ORI", "ORI: %d \n", ori);
		 p_ctg_seq("QUERY", query);
		 p_query("MATE", tmp);
		 p_query("ORIG", mate);
		 p_query("USED", get_mate(mate, seqs));
		 show_debug_msg(__func__, "OVERLAP 1: %d \n", overlapped);
		 }*/
		if (overlapped >= ol) {
			// Only if this mate overlaps with some read in the cur_pool, add it.
			// It is important because sometimes it maybe added just for coincidence.
			mate->cursor = ori ? (mate->len - overlapped - 1) : overlapped;
			pool_add(cur_pool, mate, eg->tid);
			//show_debug_msg(__func__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
			//show_debug_msg(__func__, "ORI%d \n", ori);
			//p_query("ADDED MATE", mate);
			//p_ctg_seq("TEMPLATE", template);
			//show_debug_msg(__func__, "Overlapped: %d \n", overlapped);
			//show_debug_msg(__func__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
			rm_read_from_ht(rht, mate);
		}
		if (tmp != mate)
			bwa_free_read_seq(1, tmp);
	}
	if (mate_pool->n >= N_BIG_MATE_POOL)
		g_ptr_array_free(ol_mates, TRUE);
	bwa_free_read_seq(1, template);
	destroy_reads_ht(rht);
	free_mate_pool(mate_pool);
}

