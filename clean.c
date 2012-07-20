#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include <math.h>
#include "bwase.h"
#include "peseq.h"
#include "utils.h"
#include "clean.h"

int clean_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta clean [options] \n");
	return 1;
}

double std_dev(double a[], int n) {
	if (n == 0)
		return 0.0;
	double sum = 0;
	double sq_sum = 0;
	int i = 0;
	for (i = 0; i < n; ++i) {
		sum += a[i];
		sq_sum += a[i] * a[i];
	}
	double mean = sum / n;
	double variance = sq_sum / n - mean * mean;
	return sqrt(variance);
}

double mean(double a[], int n) {
	int i = 0;
	double sum = 0;
	for (i = 0; i < n; i++) {
		sum += a[i];
	}
	return sum / n;
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

void set_kmer_index(const bwa_seq_t *read, int k, uint16_t *kmer_list) {
	int i = 0, start = 0;
	int key = 0;
	for (start = 0; start <= read->len - k; start++) {
		key = 0;
		for (i = 0; i < k; i++) {
			key *= 4;
			key = key | read->seq[i + start];
			if (read->seq[i + start] >= 4)
				return;
		}
		kmer_list[key]++;
		//		printf("Kmer: %d -> %d \n", key, kmer_list[key]);
	}
}

/* qsort int comparison function */
int cmp_kmer(const void *a, const void *b) {
	const counter *ia = (const counter *) a; // casting pointer types
	const counter *ib = (const counter *) b;
	return ib->k_freq - ia->k_freq;
}

int cmp_int(const void *a, const void *b) {
	return (*(uint16_t*) b - *(uint16_t*) a);
}

void pe_clean_core(char *fa_fn, clean_opt *opt) {
	bwa_seq_t *seqs, *part_seqs, *s, *s_pre;
	bwa_seqio_t *ks;
	int i = 0, j = 0, key = 0, tmp = 0,	zero_point = 0;
	uint32_t n_seqs = 0, n_part_seqs = 0, n_seqs_full = 0;
	int n_dup = 0, n_bad = 0, n_solid = 0, dup_thre = 0;
	uint32_t n_kmers = 0;
	char dist_name[BUFSIZE], item[BUFSIZE], solid[BUFSIZE];
	FILE *dist_file = 0, *solid_file;
	clock_t t = clock();
	uint16_t *kmer_list, *kmer_ordered, fre_upper = 0, fre_lower = 0;
	counter *k_count = 0, *k_count_pre = 0, *counter_list = 0;
	double kmer_for_read[BUFSIZE];

	ks = bwa_open_reads(opt->mode, fa_fn);
	sprintf(dist_name, "%s.kmer.dist", opt->lib_name);
	dist_file = xopen(dist_name, "w");
	sprintf(solid, "%s.solid", opt->lib_name);
	solid_file = xopen(solid, "w");

	show_debug_msg(__func__, "Loading library %s...\n", fa_fn);
	n_seqs_full = N_CHUNK_SEQS;
	seqs = (bwa_seq_t*) calloc (N_DF_MAX_SEQS, sizeof(bwa_seq_t));
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, opt->mode, 0))
			!= 0) {
		fprintf(stderr,
				"[pe_clean_core] %d sequences in library %s loaded: %.2f sec... \n",
				n_seqs + n_part_seqs, fa_fn, (float) (clock() - t) / CLOCKS_PER_SEC);
		pe_reverse_seqs(part_seqs, n_part_seqs);

		if ((n_seqs + n_part_seqs) > n_seqs_full) {
			n_seqs_full += n_part_seqs + 2;
			kroundup32(n_seqs_full);
			seqs = (bwa_seq_t*) realloc(seqs, sizeof(bwa_seq_t) * n_seqs_full);
		}
		memmove(&seqs[n_seqs], part_seqs, sizeof(bwa_seq_t) * n_part_seqs);
		free(part_seqs);
		n_seqs += n_part_seqs;
	}
	if (n_seqs < 1) {
		err_fatal(__func__,
				"No sequence in file %s, make sure the format is correct! \n",
				fa_fn);
	}
	counter_list = (counter*) calloc(n_seqs, sizeof(counter));
	n_solid = n_seqs;

	n_kmers = (1 << (opt->kmer * 2)) + 1;
	if (opt->kmer >= 16)
		n_kmers = MAX_N16;
	show_debug_msg(__func__, "# of kmers: %d\n", n_kmers);
	kmer_list = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	kmer_ordered = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	for (i = 0; i < n_kmers; i++) {
		kmer_list[i] = 0;
		kmer_ordered[i] = 0;
	}
	show_debug_msg(__func__, "Calculating k-mer frequency...\n");
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		set_kmer_index(s, opt->kmer, kmer_list);
	}
	memcpy(kmer_ordered, kmer_list, n_kmers);
//	for (i = 0; i < n_kmers; i++) {
//		kmer_ordered[i] = kmer_list[i];
//		sprintf(item, "%d\n", kmer_list[i]);
//		fputs(item, dist_file);
//	}
	show_debug_msg(__func__, "Sorting array of k-mer frequencies...\n");
	qsort(kmer_ordered, n_kmers, sizeof(uint16_t), cmp_int);
	show_debug_msg(__func__, "Getting k_freq upper and lower bounds...\n");
	for (i = 0; i < n_kmers; i++) {
		if (kmer_ordered[i] < 5) {
			zero_point = i;
			show_debug_msg(__func__, "There are %d 0-kmers.\n",
					(n_kmers - zero_point));
			break;
		}
	}
	tmp = ((zero_point) * 0.005);
	show_debug_msg(__func__, "Upper index k_freq: %d\n", tmp);
	fre_upper = kmer_ordered[tmp];
	tmp = ((zero_point) * 0.6);
	show_debug_msg(__func__, "Lower index k_freq: %d\n", tmp);
	fre_lower = kmer_ordered[tmp];
	show_debug_msg(__func__, "Lower bound k_freq: %d\n", fre_lower);
	show_debug_msg(__func__, "Upper bound k_freq: %d\n", fre_upper);
	show_debug_msg(__func__, "Ignoring reads with low frequence kmers...\n");
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		k_count = &counter_list[i];
		k_count->read_id = atoi(s->name);
		for (j = 0; j <= s->len - opt->kmer; j++) {
			key = get_key(s, j, j + opt->kmer);
			if (kmer_list[key] < fre_lower || kmer_list[key] > fre_upper) {
				k_count->k_freq = -1;
				n_bad++;
				break;
			} else {
				k_count->k_freq += kmer_list[key];
			}
		}
		if ((s->len - opt->kmer + 1) > 0)
			k_count->k_freq /= (s->len - opt->kmer + 1);
	}
	qsort(counter_list, n_solid, sizeof(counter), cmp_kmer);
	n_solid -= n_bad;
	show_debug_msg(__func__, "%d bad reads ignored, %d remained.\n", n_bad,
			n_solid);
	show_debug_msg(__func__, "Ignoring duplicate reads...\n");
	for (i = 1; i < n_solid; i++) {
		k_count = &counter_list[i];
		s = &seqs[k_count->read_id];
		k_count_pre = &counter_list[i - 1];
		j = i - 1;
		while (k_count_pre->k_freq == -1)
			k_count_pre = &counter_list[--j];
		s_pre = &seqs[k_count_pre->read_id];
		if (k_count->k_freq == k_count_pre->k_freq && same_q(s, s_pre)) {
			k_count->k_freq = -1;
			n_dup++;
		}
	}
	qsort(counter_list, n_solid, sizeof(counter), cmp_kmer);
	n_solid -= n_dup;
	show_debug_msg(__func__, "%d duplicates, %d bad, %d reads remain. \n",
			n_dup, n_bad, n_solid);
	show_debug_msg(__func__, "Saving k-mer frequencies...\n");
	for (i = 0; i < n_solid; i++) {
		k_count = &counter_list[i];
		if (k_count->k_freq > 0) {
			sprintf(item, "%d\n", k_count->read_id);
			fputs(item, solid_file);
		}
	}

	bwa_seq_close(ks);
	fclose(dist_file);
	fclose(solid_file);
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
	fprintf(stderr, "[clean_reads] Cleaning done: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return 0;
}
