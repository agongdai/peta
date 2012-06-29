#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include "bwase.h"
#include "peseq.h"
#include "utils.h"
#include "clean.h"

int clean_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta clean [options] \n");
	return 1;
}

clean_opt *init_clean_opt() {
	clean_opt *o = (clean_opt*) malloc(sizeof(clean_opt));
	o->kmer = 0;
	o->lib_name = 0;
	o->mode = (BWA_MODE_GAPE | BWA_MODE_COMPREAD);
	return o;
}

void set_kmer_index(bwa_seq_t *read, int kmer, int *kmer_list) {
	int i = 0, start = 0;
	index64 key = 0;
	for (start = 0; start <= read->len - kmer; start++) {
		key = 0;
		for (i = 0; i < kmer; i++) {
			key <<= 2;
			key = key | read->seq[i + start];
		}
		kmer_list[key]++;
	}
}

void pe_clean_core(char *fa_fn, clean_opt *opt) {
	bwa_seq_t *seqs, *s;
	bwa_seqio_t *ks;
	index64 i = 0, n_seqs = 0, n_kmers = 0;
	char dist_name[BUFSIZE], item[BUFSIZE];
	FILE *dist_file = 0;
	int *kmer_list = 0;
	clock_t t = clock();

	ks = bwa_open_reads(opt->mode, fa_fn);
	sprintf(dist_name, "%s.kmer.dist", opt->lib_name);
	dist_file = xopen(dist_name, "w");
	show_debug_msg(__func__, "Loading library %s...\n", fa_fn);
	while ((seqs = bwa_read_seq(ks, 0xa00000, &n_seqs, opt->mode, 0)) != 0) {
		fprintf(stderr,
				"[pe_clean_core] %" ID64 " sequences in library %s loaded: %.2f sec... \n",
				n_seqs, fa_fn, (float) (clock() - t) / CLOCKS_PER_SEC);
		if (n_seqs < 1) {
			err_fatal(
					__func__,
					"No sequence in file %s, make sure the format is correct! \n",
					fa_fn);
		}
		break;
	}

	n_kmers = (1 << (opt->kmer * 2)) + 1;
	show_debug_msg(__func__, "# of kmers: %" ID64 "\n", n_kmers);
	kmer_list = (int*) calloc(n_kmers, sizeof(int));
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		show_debug_msg(__func__, "Read %s length %d\n", s->name, s->len);
		set_kmer_index(s, opt->kmer, kmer_list);
	}
	for (i = 0; i < n_kmers; i++) {
		sprintf(item, "%" ID64 ": %d\n", i, kmer_list[i]);
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
