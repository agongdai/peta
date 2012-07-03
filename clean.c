#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include "bwase.h"
#include "peseq.h"
#include "utils.h"
#include "clean.h"

int clean_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta clean [options] \n");
	return 1;
}

counter *init_counter() {
	counter *mer = (counter*) malloc(sizeof(counter));
	mer->k_freq = 0;
	mer->read_id = 0;
	return mer;
}

clean_opt *init_clean_opt() {
	clean_opt *o = (clean_opt*) malloc(sizeof(clean_opt));
	o->kmer = 0;
	o->lib_name = 0;
	o->mode = (BWA_MODE_GAPE | BWA_MODE_COMPREAD);
	return o;
}

int get_key(bwa_seq_t *read, const int start, const int end) {
	int key = 0, i = 0;
	for (i = start; i < end; i++) {
		key *= 4;
		key = key | read->seq[i];
		if (read->seq[i] == 4)
			return -1;
	}
	return key;
}

void set_kmer_index(bwa_seq_t *read, int k, uint16_t *kmer_list) {
	int i = 0, start = 0;
	int key = 0;
	for (start = 0; start <= read->len - k; start++) {
		key = 0;
		for (i = 0; i < k; i++) {
			key *= 4;
			key = key | read->seq[i + start];
			if (read->seq[i + start] == 4)
				return;
		}
		kmer_list[key]++;
		//		printf("Kmer: %d -> %d \n", key, kmer_list[key]);
	}
}

/* qsort int comparison function */
int int_cmp(const void *a, const void *b) {
	const counter *ia = (const counter *) a; // casting pointer types
	const counter *ib = (const counter *) b;
	return ib->k_freq - ia->k_freq;
}

void pe_clean_core(char *fa_fn, clean_opt *opt) {
	bwa_seq_t *seqs, *s;
	bwa_seqio_t *ks;
	int i = 0, j = 0, n_seqs = 0, n_kmers = 0, key = 0;
	char dist_name[BUFSIZE], item[BUFSIZE];
	FILE *dist_file = 0;
	clock_t t = clock();
	uint16_t *kmer_list;
	uint16_t mer = 0;
	counter *k_count = 0, *counter_list = 0;

	ks = bwa_open_reads(opt->mode, fa_fn);
	sprintf(dist_name, "%s.kmer.dist", opt->lib_name);
	dist_file = xopen(dist_name, "w");
	show_debug_msg(__func__, "Loading library %s...\n", fa_fn);
	while ((seqs = bwa_read_seq(ks, 0xa00000, &n_seqs, opt->mode, 0)) != 0) {
		fprintf(
				stderr,
				"[pe_clean_core] %d sequences in library %s loaded: %.2f sec... \n",
				n_seqs, fa_fn, (float) (clock() - t) / CLOCKS_PER_SEC);
		pe_reverse_seqs(seqs, n_seqs);
		counter_list = (counter*) calloc(n_seqs, sizeof(counter));
		if (n_seqs < 1) {
			err_fatal(
					__func__,
					"No sequence in file %s, make sure the format is correct! \n",
					fa_fn);
		}
		break;
	}

	n_kmers = (1 << (opt->kmer * 2)) + 1;
	show_debug_msg(__func__, "# of kmers: %d\n", n_kmers);
	kmer_list = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	for (i = 0; i < n_kmers; i++) {
		kmer_list[i] = 0;
	}
	show_debug_msg(__func__, "Calculating k-mer frequency...\n");
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		//		p_query(__func__, s);
		set_kmer_index(s, opt->kmer, kmer_list);
	}
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		k_count = &counter_list[i];
		k_count->read_id = s->name;
		for (j = 0; j <= s->len - opt->kmer; j++) {
			key = get_key(s, j, j + opt->kmer);
			k_count->k_freq += kmer_list[key];
		}
	}
	show_debug_msg(__func__, "Sorting array of k-mer frequencies...\n");
	qsort(counter_list, n_seqs, sizeof(counter), int_cmp);
	for (i = n_kmers - 1; i > n_kmers - 10; i--) {
		mer = kmer_list[i];
		printf("%d coverage: %d \n", i, mer);
	}
	show_debug_msg(__func__, "Saving k-mer frequencies...\n");
	for (i = 0; i < n_seqs; i++) {
		k_count = &counter_list[i];
		sprintf(item, "%s\n", k_count->read_id);
		fputs(item, dist_file);
	}

	bwa_seq_close(ks);
	fclose(dist_file);
}

int clean_reads(int argc, char *argv[]) {
	int c;
	clock_t t;
	t = clock();
	clean_opt *opt = init_clean_opt();

	while ((c = getopt(argc, argv, "k:l:")) >= 0) {
		switch (c) {
		case 'k':
			opt->kmer = atoi(optarg);
			break;
		case 'l':
			opt->lib_name = optarg;
			break;
		default:
			return 1;
		}
	}

	if (optind + 1 > argc) {
		return clean_usage();
	}

	pe_clean_core(argv[optind], opt);
	fprintf(stderr, "[clean_reads] Cleaning done: %.2f sec\n", (float) (clock()
			- t) / CLOCKS_PER_SEC);
	return 0;
}
