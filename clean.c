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
#include "pehash.h"
#include "rnaseq.h"

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
	mer->k_freq = 0.0;
	mer->read_id = 0;
	mer->k_sd = 0.0;
	mer->k_mean = 0.0;
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

int has_low_kmer(bwa_seq_t *read, uint16_t *kmer_list, clean_opt *opt) {
	int j = 0, key = 0;
	for (j = 0; j <= read->len - opt->kmer; j++) {
		key = get_key(read, j, j + opt->kmer);
		if (kmer_list[key] <= LOW_KMER) {
			return 1;
		}
	}
	return 0;
}

int kmer_uneven(bwa_seq_t *read, counter *counter, uint16_t *kmer_list,
		clean_opt *opt, double range) {
	int i = 0, key = 0, j = 0;
	double lower_bound = 0, upper_bound = 0;
	// For every base in the read, count how many kmers covering it.
	double *base_counter = (double*) calloc(read->len + 1, sizeof(double));
	for (i = 0; i <= read->len - opt->kmer; i++) {
		key = get_key(read, j, j + opt->kmer);
		for (j = i; j < opt->kmer + i; j++) {
			base_counter[j] += kmer_list[key];
		}
	}
	show_debug_msg(__func__, "Counted \n");
	counter->k_mean = mean(base_counter, read->len);
	show_debug_msg(__func__, "Mean: %.2f \n", counter->k_mean);
	counter->k_sd = std_dev(base_counter, read->len);
	show_debug_msg(__func__, "SD: %.2f \n", counter->k_sd);
	lower_bound = counter->k_mean - counter->k_sd * range;
	upper_bound = counter->k_mean + counter->k_sd * range;
	show_debug_msg(__func__, "Range: [%.2f, %.2f] \n", lower_bound, upper_bound);
	for (i = 0; i < read->len; i++) {
		// If any base in the read has lower or higher coverage then the mean, the read is uneven
		if (base_counter[i] > upper_bound || base_counter[i] < lower_bound) {
			return 1;
		}
	}
	return 0;
}

void pe_clean_core(char *fa_fn, clean_opt *opt) {
	bwa_seq_t *seqs, *s, *s_pre;
	int i = 0, j = 0, key = 0, tmp = 0, zero_point = 0;
	int n_dup = 0, n_bad = 0, n_solid = 0, dup_thre = 0, n_uneven = 0;
	uint32_t n_kmers = 0;
	char item[BUFSIZE], solid[BUFSIZE];
	FILE *dist_file = 0, *solid_file;
	clock_t t = clock();
	uint16_t *kmer_list, *kmer_ordered, fre_upper = 0, fre_lower = 0;
	counter *k_count = 0, *k_count_pre = 0, *counter_list = 0;
	double kmer_for_read[BUFSIZE];
	hash_table *ht = NULL;

	sprintf(solid, "%s.solid", opt->lib_name);
	solid_file = xopen(solid, "w");

	show_debug_msg(__func__, "Loading library %s...\n", fa_fn);
	ht = pe_load_hash(fa_fn);
	seqs = ht->seqs;
	// Each counter corresponds to a read.
	// It contains the kmer frequencies.
	counter_list = (counter*) calloc(ht->n_seqs, sizeof(counter));
	n_solid = ht->n_seqs;

	n_kmers = (1 << (opt->kmer * 2)) + 1;
	if (opt->kmer >= 16)
		n_kmers = MAX_N16;
	show_debug_msg(__func__, "kmer length: %d; # of kmers: %d\n", opt->kmer,
			n_kmers);
	kmer_list = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	kmer_ordered = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	for (i = 0; i < n_kmers; i++) {
		kmer_list[i] = 0;
		kmer_ordered[i] = 0;
	}

	// For every kmer in the read, count it and store it in kmer_list
	show_debug_msg(__func__, "Calculating k-mer frequency: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 0; i < ht->n_seqs; i++) {
		s = &seqs[i];
		set_kmer_index(s, opt->kmer, kmer_list);
	}
	// kmer_ordered is a clone of kmer_list, sorted.
	memcpy(kmer_ordered, kmer_list, n_kmers);

	show_debug_msg(__func__,
			"Sorting array of k-mer frequencies: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	qsort(kmer_ordered, n_kmers, sizeof(uint16_t), cmp_int);
	show_debug_msg(
			__func__,
			"Ignoring reads having kmers with low kmer frequencies: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 0; i < ht->n_seqs; i++) {
		s = &seqs[i];
		k_count = &counter_list[i];
		k_count->read_id = atoi(s->name);
		if (has_low_kmer(s, kmer_list, opt)) {
			k_count->k_freq = -1;
			n_bad++;
		} else {
			if (kmer_uneven(read, k_count, kmer_list, opt, UNEVEN_THRE)) {
				k_count->k_freq = -1;
				n_uneven++;
			}
		}
	}
	qsort(counter_list, n_solid, sizeof(counter), cmp_kmer);
	n_solid -= n_bad;
	n_solid -= n_uneven;
	show_debug_msg(__func__, "%d bad reads, % uneven reads, %d remained.\n",
			n_bad, n_uneven, n_solid);
	show_debug_msg(__func__, "Ignoring duplicate reads: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 1; i < n_solid; i++) {
		k_count = &counter_list[i];
		s = &seqs[k_count->read_id];
		k_count_pre = &counter_list[i - 1];
		j = i - 1;
		while (k_count_pre->k_freq == -1)
			k_count_pre = &counter_list[--j];
		s_pre = &seqs[k_count_pre->read_id];
		// If two reads has the same kmer frequencies and their sequences are the same, remove it.
		if (k_count->k_freq == k_count_pre->k_freq && same_q(s, s_pre)) {
			k_count->k_freq = -1;
			n_dup++;
		}
	}
	qsort(counter_list, n_solid, sizeof(counter), cmp_kmer);
	n_solid -= n_dup;
	show_debug_msg(__func__, "%d duplicates, %d bad, %d reads remain. \n",
			n_dup, n_bad, n_solid);
	show_debug_msg(__func__, "Saving k-mer frequencies: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 0; i < n_solid; i++) {
		k_count = &counter_list[i];
		if (k_count->k_freq > 0) {
			sprintf(item, "%d\n", k_count->read_id);
			fputs(item, solid_file);
		}
	}
	show_debug_msg(__func__, "%d solid reads remained: %.2f.\n", n_solid,
			(float) (clock() - t) / CLOCKS_PER_SEC);
	fclose(solid_file);
}

void test_clean() {
	double numbers[] = { 2.0, 4.0, 4.0, 4.0, 5, 5.0, 7.0, 9.0 };
	printf("Mean: %.3f \n", mean(numbers, 8));
	printf("SD: %.3f \n", std_dev(numbers, 8));
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
	//	test_clean();
	fprintf(stderr, "[clean_reads] Cleaning done: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return 0;
}
