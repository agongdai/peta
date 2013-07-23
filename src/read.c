#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include "k_hash.h"
#include "bwtaln.h"
#include "peseq.h"
#include "rnaseq.h"
#include "read.h"

int group_usage() {
	printf("Usage of group: \n");
	return 1;
}

gint cmp_hits_by_count(gpointer a, gpointer b) {
	GPtrArray *hit_a = *((GPtrArray**) a);
	GPtrArray *hit_b = *((GPtrArray**) b);
	return (hit_b->len - hit_a->len);
}

read_hash *new_rh(bwa_seq_t *seqs, uint64_t n_seqs) {
	GPtrArray *hits = NULL;
	int i = 0;
	read_hash *rh = (read_hash*) malloc(sizeof(read_hash));
	rh->seqs = seqs;
	rh->n_seqs = n_seqs;
	rh->read_groups = g_ptr_array_sized_new(n_seqs);
	for (i = 0; i < n_seqs; i++) {
		hits = g_ptr_array_sized_new(0);
		g_ptr_array_add(rh->read_groups, hits);
	}
	return rh;
}

void destroy_rh(read_hash *rh) {
	int i = 0;
	GPtrArray *hits = NULL;
	if (!rh)
		return;
	show_msg(__func__, "Destroying read hash ...\n");
	if (rh->read_groups) {
		for (i = 0; i < rh->read_groups->len; i++) {
			hits = (GPtrArray*) g_ptr_array_index(rh->read_groups, i);
			if (hits) {
				g_ptr_array_free(hits, TRUE);
				hits = NULL;
			}
		}
		g_ptr_array_free(rh->read_groups, TRUE);
	}
	if (rh->seqs && rh->n_seqs > 0) {
		bwa_free_read_seq(rh->n_seqs, rh->seqs);
	}
	free(rh);
}

void p_read_hash(read_hash *rh) {
	int i = 0, j = 0;
	bwa_seq_t *query = NULL, *r = NULL;
	GPtrArray *hits = NULL;
	show_msg(__func__, "Printing groups ...\n");
	for (i = 0; i < rh->read_groups->len; i++) {
		query = &rh->seqs[i];
		p_shift_query(query, query->len);
		hits = (GPtrArray*) g_ptr_array_index(rh->read_groups, i);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			p_shift_query(r, r->len);
		}
		show_debug_msg(__func__, "----------------------------------------\n");
		if (i > 20)
			break;
	}
}

/**
 * Save the read groups as binary file
 */
void dump_read_hash(char *fa, read_hash *rh) {
	char dump_fn[BUFSIZ];
	FILE *group = NULL;
	uint64_t id = 0, *hit_ids = NULL;
	int i = 0, j = 0;
	bwa_seq_t *r = NULL;
	GPtrArray *read_groups = rh->read_groups, *hits = NULL;
	sprintf(dump_fn, "%s.group", fa);
	group = xopen(dump_fn, "w");

	show_msg(__func__, "Dumping read groups to file %s ...\n", dump_fn);
	for (i = 0; i < read_groups->len; i++) {
		hits = (GPtrArray*) g_ptr_array_index(read_groups, i);
		hit_ids = (uint64_t*) calloc(hits->len + 1, sizeof(uint64_t));
		hit_ids[0] = hits->len;
		//r = &rh->seqs[i];
		//p_query(__func__, r);
		for (j = 0; j < hits->len; j++) {
			r = (bwa_seq_t*) g_ptr_array_index(hits, j);
			//p_query(__func__, r);
			hit_ids[j + 1] = atol(r->name);
		}
		fwrite(hit_ids, sizeof(uint64_t), hits->len + 1, group);
	}
	show_msg(__func__, "Read groups dumped to file %s \n", dump_fn);
	fclose(group);
}

read_hash *load_read_hash(char *fa) {
	uint64_t n_seqs = 0;
	bwa_seq_t *seqs = load_reads(fa, &n_seqs);
	bwa_seq_t *r = NULL;
	read_hash *rh = new_rh(seqs, n_seqs);
	GPtrArray *hits = NULL;
	uint64_t *hit_ids = NULL, id = 0, n_hits = 0;
	char dump_fn[BUFSIZ];
	FILE *group = NULL;
	sprintf(dump_fn, "%s.group", fa);
	group = xopen(dump_fn, "r");

	show_msg(__func__, "Loading read groups from %s ...\n", dump_fn);
	int i = 0, j = 0;
	rh->seqs = seqs;
	rh->n_seqs = n_seqs;
	for (i = 0; i < rh->read_groups->len; i++) {
		r = &seqs[i];
		hits = (GPtrArray*) g_ptr_array_index(rh->read_groups, i);
		fread(&n_hits, sizeof(uint64_t), 1, group);
		//show_debug_msg(__func__, "%d hits for read %s \n", n_hits, r->name);
		hit_ids = (uint64_t*) calloc(n_hits, sizeof(uint64_t));
		fread(hit_ids, sizeof(uint64_t), n_hits, group);
		for (j = 0; j < n_hits; j++) {
			id = hit_ids[j];
			r = &seqs[id];
			g_ptr_array_add(hits, r);
		}
	}
	return rh;
}

/**
 * The thread to align one read to the RNA-seq library
 */
gint align_read_thread(gpointer r, gpointer para) {
	bwa_seq_t *query = (bwa_seq_t*) r;
	align_para *p = (align_para*) para;
	hash_table *ht = p->ht;
	GPtrArray *hits = NULL;
	GPtrArray *read_groups = p->rh->read_groups;
	if (atoi(query->name) < 0 || atoi(query->name) >= ht->n_seqs)
		return 1;
	hits = (GPtrArray*) g_ptr_array_index(read_groups, atoi(query->name));
	find_both_fr_reads(ht, query, hits, N_MISMATCHES);
	return 0;
}

/**
 * Start n_threads threads to align reads
 */
read_hash *group_reads(hash_table *ht, const int n_threads) {
	GThreadPool *thread_pool = NULL;
	int i = 0;
	bwa_seq_t *r = NULL;
	read_hash *rh = new_rh(ht->seqs, ht->n_seqs);
	GPtrArray *read_groups = rh->read_groups;
	align_para *para = (align_para*) malloc(sizeof(align_para));

	para->ht = ht;
	para->rh = rh;

	show_msg(__func__, "Aligning %d reads on %d threads ...\n", ht->n_seqs, n_threads);
	thread_pool = g_thread_pool_new((GFunc) align_read_thread, para, n_threads,
			TRUE, NULL);
	for (i = 0; i < ht->n_seqs; i++) {
		r = (bwa_seq_t*) &ht->seqs[i];
		g_thread_pool_push(thread_pool, (gpointer) r, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	free(para);

	return rh;
}

void test_group(char *fa, int n_threads) {
	char *fa_copy = strdup(fa);
	char *fa_copy2 = strdup(fa);
	hash_table *ht = load_k_hash(fa);
	read_hash *rh = group_reads(ht, n_threads);
	dump_read_hash(fa_copy, rh);

	p_read_hash(rh);
	destroy_rh(rh);
	rh = load_read_hash(fa_copy2);
	show_debug_msg(__func__, "TAG: Here are hits loaded\n");
	p_read_hash(rh);
	destroy_rh(rh);
	free(fa_copy);
}

int group_reads_core(char *fa, int n_threads) {
	char *fa_copy = strdup(fa);
	char *fa_copy2 = strdup(fa);
	hash_table *ht = load_k_hash(fa_copy);
	read_hash *rh = group_reads(ht, n_threads);
	dump_read_hash(fa_copy2, rh);
	free(fa_copy);
	free(fa_copy2);
}

int group_main(int argc, char *argv[]) {
	int c;
	clock_t t;
	t = clock();
	int n_threads = 1;
	char *fa = NULL;

	while ((c = getopt(argc, argv, "t:")) >= 0) {
		switch (c) {
		case 't':
			n_threads = atoi(optarg);
			break;
		default:
			return 1;
		}
	}

	if (optind + 1 > argc) {
		return group_usage();
	}

	//test_group(argv[optind], n_threads);
	fa = strdup(argv[optind]);
	group_reads_core(fa, n_threads);
	free(fa);
	fprintf(stderr, "[group_main] Grouping done: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return 0;
}
