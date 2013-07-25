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

read_hash *new_rh(uint64_t n_seqs) {
	GPtrArray *hits = NULL;
	int i = 0;
	read_hash *rh = (read_hash*) malloc(sizeof(read_hash));
	rh->n_seqs = n_seqs;
	rh->similar_reads_count = (index64*) calloc(n_seqs, sizeof(index64));
	return rh;
}

void p_read_hash(read_hash *rh) {

}

void destroy_rh(read_hash *rh) {
	int i = 0;
	GPtrArray *hits = NULL;
	if (!rh)
		return;
	if (rh->similar_reads_count) {
		free(rh->similar_reads_count);
	}
	free(rh);
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
	index64 *similar_reads_count = rh->similar_reads_count, *hits = NULL;
	sprintf(dump_fn, "%s.group", fa);
	group = xopen(dump_fn, "w");

	show_msg(__func__, "Dumping read groups to file %s ...\n", dump_fn);
	fwrite(&rh->n_seqs, sizeof(index64), 1, group);
	fwrite(rh->similar_reads_count, sizeof(index64), rh->n_seqs, group);
	show_msg(__func__, "Similar read counts dumped to file %s. \n", dump_fn);
	fclose(group);
}

read_hash *load_read_hash(char *fa) {
	uint64_t n_seqs = 0;
	char dump_fn[BUFSIZ];
	FILE *group = NULL;
	sprintf(dump_fn, "%s.group", fa);
	group = xopen(dump_fn, "r");

	show_msg(__func__, "Loading read groups from %s ...\n", dump_fn);
	fread(&n_seqs, sizeof(index64), 1, group);
	read_hash *rh = new_rh(n_seqs);
	fread(rh->similar_reads_count, sizeof(uint64_t), n_seqs, group);
	fclose(group);
	return rh;
}

/**
 * The thread to align one read to the RNA-seq library
 */
gint align_read_thread(gpointer r, gpointer para) {
	bwa_seq_t *query = (bwa_seq_t*) r;
	align_para *p = (align_para*) para;
	hash_table *ht = p->ht;
	GPtrArray *hits = g_ptr_array_sized_new(0);
	index64 i = 0;
	bwa_seq_t *tmp = NULL;
	//p_query(__func__, query);
	//show_debug_msg(__func__, "# of n_seqs: %d\n", ht->n_seqs);
	index64 *similar_reads_count = p->rh->similar_reads_count;
	index64 read_id = 0;
	read_id = atoll(query->name);
	if (read_id < 0 || read_id >= ht->n_seqs || query->pos != -1)
		return 1;
	find_both_fr_full_reads(ht, query, hits, N_MISMATCHES);
	// Mark 'pos' as not -1, to save time
	for (i = 0; i < hits->len; i++) {
		tmp = (bwa_seq_t*) g_ptr_array_index(hits, i);
		tmp->pos = hits->len;
		similar_reads_count[atoll(tmp->name)] = hits->len;
	}
	//show_debug_msg(__func__, "Hits: %d\n", hits->len);
	//show_debug_msg(__func__, "Query id: %d\n", read_id);
	similar_reads_count[read_id] = hits->len;
	g_ptr_array_free(hits, TRUE);
	return 0;
}

/**
 * Start n_threads threads to align reads
 */
read_hash *group_reads(hash_table *ht, const int n_threads) {
	GThreadPool *thread_pool = NULL;
	int i = 0;
	bwa_seq_t *r = NULL;
	read_hash *rh = new_rh(ht->n_seqs);
	align_para *para = (align_para*) malloc(sizeof(align_para));

	para->ht = ht;
	para->rh = rh;

	show_msg(__func__, "Aligning %d reads on %d threads ...\n", ht->n_seqs,
			n_threads);
	thread_pool = g_thread_pool_new((GFunc) align_read_thread, para, n_threads,
			TRUE, NULL);
	for (i = 0; i < ht->n_seqs; i++) {
		r = (bwa_seq_t*) &ht->seqs[i];
		//align_read_thread((gpointer)r, (gpointer)para);
		g_thread_pool_push(thread_pool, (gpointer) r, NULL);
	}
	g_thread_pool_free(thread_pool, 0, 1);
	free(para);

	return rh;
}

int group_reads_core(char *fa, int n_threads) {
	char *fa_copy = strdup(fa);
	char *fa_copy2 = strdup(fa);
	hash_table *ht = load_k_hash(fa_copy);
	read_hash *rh = group_reads(ht, n_threads);
	dump_read_hash(fa_copy2, rh);
	destroy_ht(ht);
	free(fa_copy);
	free(fa_copy2);
}

void test_group(char *fa, int n_threads) {

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

	if (!g_thread_supported())
		g_thread_init(NULL);
	//test_group(argv[optind], n_threads);
	fa = strdup(argv[optind]);
	group_reads_core(fa, n_threads);
	free(fa);
	fprintf(stderr, "[group_main] Grouping done: %.2f sec\n", (float) (clock()
			- t) / CLOCKS_PER_SEC);
	return 0;
}
