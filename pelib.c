#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <math.h>
#include "rand.h"
#include "ass.h"
#include "bwase.h"
#include "utils.h"
#include "peseq.h"
#include "pool.h"
#include "pechar.h"
#include "pealn.h"
#include "roadmap.h"
#include "edge.h"
#include "edgelist.h"
#include "readrm.h"
#include "pelib.h"

void pe_lib_help() {
	show_msg(__func__, "--------------------------------------------------");
	show_msg(__func__, "peta pelib <library>");
	show_msg(__func__, "--------------------------------------------------");
}

void maintain_pool(alignarray *aligns, const hash_table *ht, pool *cur_pool,
		pool *mate_pool, bwa_seq_t *contig, bwa_seq_t *query, int *next) {
	alg *a = NULL;
	int i = 0, index = 0;
	bwa_seq_t *s = NULL, *mate = NULL, *tmp = NULL, *seqs = ht->seqs;
	for (i = 0; i < aligns->len; i++) {
		a = g_ptr_array_index(aligns, i);
		index = a->r_id;
		if (index >= ht->n_seqs)
			continue;
		s = &seqs[index];
		s->rev_com = a->rev_comp;
		mate = get_mate(s, seqs);
		mate->rev_com = s->rev_com;
		if (s->is_in_c_pool || (a->pos + query->len - 1) > s->len || s->used)
			continue;
		s->cursor = s->rev_com ? (s->len - a->pos) : (a->pos + query->len);
		if (s->cursor >= s->len || s->cursor < 0) {
			s->cursor = 0;
			continue;
		}
		pool_add(cur_pool, s);
		if (is_right_mate(mate->name) && !mate->used && !mate->is_in_c_pool)
			pool_add(mate_pool, mate);
	}
	// Add the mate reads which overlap with the tail into the current pool
	for (i = 0; i < mate_pool->n; i++) {
		mate = g_ptr_array_index(mate_pool->reads, i);
		tmp = mate;
		if (mate->rev_com)
			tmp = new_mem_rev_seq(mate, mate->len, 0);
		if (seq_ol(contig, tmp, mate->len / 4, MISMATCHES)) {
			mate->cursor = mate->len / 4;
			pool_add(cur_pool, mate);
		}
		if (mate->rev_com)
			bwa_free_read_seq(1, tmp);
	}
	for (i = 0; i < cur_pool->n; i++) {
		s = g_ptr_array_index(cur_pool->reads, i);
		if (s->rev_com)
			check_c(next, s->rseq[s->cursor]);
		else
			check_c(next, s->seq[s->cursor]);
	}
}

edge *est_ins_size(const hash_table *ht, const bwa_seq_t *s, double *sizes,
		int *n_pairs) {
	bwa_seq_t *query = NULL;
	pool *cur_pool = NULL, *mate_pool = NULL;
	int c = 0;
	int *next = (int*) calloc(5, sizeof(int));
	alignarray *aligns = NULL;
	edge *eg = NULL;

	eg = new_eg();
	eg->contig = new_seq(s, s->len, 0);
	query = new_seq(s, s->len / 2, 0);
	cur_pool = new_pool();
	mate_pool = new_pool();
	while (1) {
		p_query(__func__, query);
		reset_c(next, NULL); // Reset the counter
		pe_aln_query(query, query->seq, ht, MISMATCHES, query->len, 0, aligns);
		pe_aln_query(query, query->rseq, ht, MISMATCHES, query->len, 1, aligns);
		maintain_pool(aligns, ht, cur_pool, mate_pool, eg->contig, query, next);
		reset_alg(aligns);
		p_pool("Current Pool", cur_pool, next);
		p_pool("Mate Pool", mate_pool, next);
		c = get_abs_most(next, STRICT_PERC);
		if (c == -1 && cur_pool->n <= 0)
			break;
		forward(cur_pool, c, eg, 0);
		ext_con(eg->contig, c, 0);
		p_ctg_seq("Contig now", eg->contig);
		eg->len = eg->contig->len;
		ext_que(query, c, 0);
	}
	free_alg(aligns);
	free(next);
	free_pool(cur_pool);
	free_pool(mate_pool);
	bwa_free_read_seq(1, query);
	return eg;
}

void pe_lib_core(char *lib_file, char *solid_file) {
	hash_table *ht = NULL;
	bwa_seq_t *query = NULL, *seqs = NULL;
	FILE *solid = xopen(solid_file, "r");
	char line[80];
	int index = 0, ol = 0, n_pairs = 0;
	double *sizes = NULL;
	edge *eg = NULL;
	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);

	sizes = (double*) calloc(ht->n_seqs / 2, sizeof(int));
	ht = pe_load_hash(lib_file);
	seqs = &ht->seqs;
	ol = seqs->len / 2; // Read length
	while (fgets(line, 80, solid) != NULL) {
		index = atoi(line);
		query = &ht->seqs[index];
		eg = est_ins_size(ht, query, sizes, &n_pairs);
		destroy_eg(eg);
		break;
	}
	free(sizes);
	fclose(solid);
	destroy_ht(ht);
}

int pe_lib(int argc, char *argv[]) {
	clock_t t = clock();
	pe_lib_core(argv[1], argv[2]);
	show_msg(__func__, "Done: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return 0;
}
