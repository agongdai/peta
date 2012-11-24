/*
 * pet.c
 *
 *  Created on: 23-Jun-2011
 *      Author: carl
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bwase.h"
#include "pool.h"
#include "peseq.h"
#include "bwtaln.h"
#include "pet.h"
#include "pechar.h"
#include "pool.h"
#include "pealn.h"
#include "pehash.h"

pool *pet_cvg(const char *pet_fn, const ass_opt *opt) {
	bwa_seq_t *pets, *query, *p, *p2;
	int i = 0, j = 0, k = 0;
	index64 mate_i = 0;
	pool *good_pets = new_pool(), *repeat_pets = new_pool();
	alignarray *align, *align_2;
	alg *a;
	hash_table *ht;

	ht = pe_load_hash(pet_fn);
	pets = ht->seqs;

	fprintf(stderr, "[pe_cvg] Converging RNA-PETs... \n");
	// for (i = n_pets - 1; i >= 0; i -= 2) {
	for (i = 0; i < ht->n_seqs; i += 2) {
		p = &pets[i];
		p2 = &pets[i + 1];
		if (binary_exists(repeat_pets->reads, p) || binary_exists(good_pets->reads, p))
			continue;

		for (k = p->len - opt->ol; k >= 0; k--) {
			query = new_seq(p, opt->ol, k);
			// p_query(query);
			pe_aln_query(query, query->seq, ht, opt->nm + 2, opt->ol, 0, align);
			pool_sort_ins(good_pets, p);
			// p_align(align);

			query = new_seq(p2, opt->ol, k);
			// p_query(query);
			pe_aln_query(query, query->seq, ht, opt->nm + 2, opt->ol, 0, align_2);
			pool_sort_ins(good_pets, p2);
			// p_align(align_2);
			for (j = 0; j < align->len; j++) {
				a = g_ptr_array_index(align, j);
				// The aligned seq is the query itself
				if (a->r_id == atoll(p->name))
					continue;
				mate_i = get_mate_index(a->r_id);
				// If the right mate is also aligned
				if (!aligned(align_2, mate_i))
					continue;
				pool_sort_ins(repeat_pets, &pets[a->r_id]);
				pool_sort_ins(repeat_pets, &pets[mate_i]);
			}
		}
		// p_pool("Good Pets: ", good_pets);
		// p_pool("Repeat Pets: ", repeat_pets);
	}

	fprintf(stderr, "[pet_cvg] Converged to %zd RNA-PETs... \n", (good_pets->n));
	fprintf(stderr, "[pet_cvg] ------------------------------ \n");
	//	p_pool("Good Pets: ", good_pets);
	return good_pets;
}
