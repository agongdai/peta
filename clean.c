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
#include "pealn.h"

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
	mer->checked = 0;
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

int pick_within_range(bwa_seq_t *read, counter *counter, uint16_t *kmer_list,
		clean_opt *opt, double range) {
	int i = 0, j = 0, key = 0;
	double lower_bound = 0, upper_bound = 0;
	double *base_counter = (double*) calloc(read->len + 1, sizeof(double));
	lower_bound = counter->k_mean - counter->k_sd * range;
	upper_bound = counter->k_mean + counter->k_sd * range;
	for (i = 0; i <= read->len - opt->kmer; i++) {
		key = get_key(read, i, i + opt->kmer);
		for (j = i; j < opt->kmer + i; j++) {
			base_counter[j] += kmer_list[key];
		}
	}
	for (i = opt->kmer; i < read->len - opt->kmer; i++) {
		// If any base in the read has lower or higher coverage then the mean, the read is uneven
		if (base_counter[i] > upper_bound || base_counter[i] < lower_bound) {
			free(base_counter);
			return 0;
		}
	}
	free(base_counter);
	return 1;
}

void set_k_freq(bwa_seq_t *read, counter *k_count, uint16_t *kmer_list,
		clean_opt *opt) {
	int j = 0, key = 0, i = 0, counted_len = read->len - 2 * opt->kmer;
	double *base_counter = (double*) calloc(counted_len + 1, sizeof(double));
	for (j = 0; j <= read->len - opt->kmer; j++) {
		key = get_key(read, j, j + opt->kmer);
		k_count->k_freq += kmer_list[key];
	}
	for (i = opt->kmer; i <= counted_len; i++) {
		key = get_key(read, i, i + opt->kmer);
		for (j = i; j < opt->kmer + i; j++) {
			base_counter[j - opt->kmer] += kmer_list[key];
		}
	}
	k_count->k_mean = mean(base_counter, counted_len);
	k_count->k_sd = std_dev(base_counter, counted_len);
	free(base_counter);
}

void mark_duplicates(bwa_seq_t *s, hash_table *ht) {
	int j = 0, significant_len = s->len - 4;
	bwa_seq_t *seqs, *part_s = NULL, *dup_seq = NULL;
	alignarray *aligns;
	alg *a;

	aligns = g_ptr_array_sized_new(N_DEFAULT_ALIGNS);
	seqs = ht->seqs;
	part_s = new_seq(s, significant_len, 0);
	pe_aln_query(part_s, part_s->seq, ht, 0, significant_len, 0, aligns);
	pe_aln_query(part_s, part_s->seq, ht, 0, significant_len, 0, aligns);
	for (j = 0; j < aligns->len; j++) {
		a = g_ptr_array_index(aligns, j);
		dup_seq = &seqs[a->r_id];
		dup_seq->used = 1;
	}
	free_alg(aligns);
	bwa_free_read_seq(1, part_s);
}

void pe_clean_core(char *fa_fn, clean_opt *opt) {
	bwa_seq_t *seqs, *s;
	int i = 0, j = 0;
	int n_dup = 0, n_bad = 0, n_solid = 0, n_rep = 0;
	uint32_t n_kmers = 0, n_seqs = 0;
	char item[BUFSIZE], solid[BUFSIZE];
	FILE *solid_file;
	clock_t t = clock();
	uint16_t *kmer_list;
	counter *k_count = NULL, *counter_list = NULL, *sorted_counters = NULL;
	GPtrArray *solid_reads = g_ptr_array_sized_new(16384);
	hash_table *ht = NULL;

	sprintf(solid, "%s.solid", opt->lib_name);
	solid_file = xopen(solid, "w");

	show_debug_msg(__func__, "Loading library %s...\n", fa_fn);
	ht = pe_load_hash(fa_fn);
	seqs = ht->seqs;
	n_seqs = ht->n_seqs;
	// Each counter corresponds to a read.
	// It contains the kmer frequencies.
	counter_list = (counter*) calloc(n_seqs, sizeof(counter));
	sorted_counters = (counter*) calloc(n_seqs, sizeof(counter));

	n_kmers = (1 << (opt->kmer * 2)) + 1;
	if (opt->kmer >= 16)
		n_kmers = MAX_N16;
	show_debug_msg(__func__, "kmer length: %d; # of kmers: %d\n", opt->kmer,
			n_kmers);
	kmer_list = (uint16_t*) calloc(n_kmers, sizeof(uint16_t));
	for (i = 0; i < n_kmers; i++) {
		kmer_list[i] = 0;
	}

	// For every kmer in the read, count it and store it in kmer_list
	show_debug_msg(__func__, "Calculating k-mer frequency: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		set_kmer_index(s, opt->kmer, kmer_list);
	}

	// Calculate the mean and standard deviation of k_freq
	show_debug_msg(__func__,
			"Counting the k-mer frequencies of reads: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		k_count = &counter_list[i];
		k_count->read_id = atoi(s->name);
		k_count->checked = 0;
		set_k_freq(s, k_count, kmer_list, opt);
	}

	show_debug_msg(__func__, "Removing repetitive and low k-mer frequency reads: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);

	// Remove repetitive reads and those reads having low frequency kmers.
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		k_count = &counter_list[i];
		if (has_rep_pattern(s)) {
			n_rep++;
			k_count->checked = 1;
		} else {
			if (has_low_kmer(s, kmer_list, opt)) {
				k_count->checked = 2;
				n_bad++;
			}
		}
	}
	show_debug_msg(__func__, "Removed %d repetitive reads, %d low k-mer reads. \n",
				n_rep, n_bad);

	// Sort the counters
	show_debug_msg(__func__, "Sorting reads by k_freq: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	memcpy(sorted_counters, counter_list, sizeof(counter) * n_seqs);
	qsort(sorted_counters, n_seqs, sizeof(counter), cmp_kmer);

	show_debug_msg(__func__,
			"Getting solid reads (duplicates removed): %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (j = 1; j <= TRIAL_TIME; j++) {
		show_debug_msg(
				__func__,
				"Round %d out of %d. %d duplicate reads; %d solid reads: %.2f sec...\n",
				j, TRIAL_TIME, n_dup, n_solid,
				(float) (clock() - t) / CLOCKS_PER_SEC);
		for (i = 0; i < n_seqs; i++) {
			k_count = &sorted_counters[i];
			s = &seqs[k_count->read_id];
			if (k_count->checked)
				continue;
			if (s->used) {
				n_dup++;
				k_count->checked = 3;
			}
			if (pick_within_range(s, k_count, kmer_list, opt, UNEVEN_THRE * j)) {
				n_solid++;
				k_count->checked = 4;
				g_ptr_array_add(solid_reads, s);
				mark_duplicates(s, ht);
			}
		}
	}

	show_debug_msg(__func__, "%d duplicate reads; %d solid reads.\n", n_dup,
			n_solid);

	show_debug_msg(__func__, "Saving k-mer frequencies: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	for (i = 0; i < solid_reads->len; i++) {
		s = g_ptr_array_index(solid_reads, i);
		sprintf(item, "%s\n", s->name);
		fputs(item, solid_file);
	}

	show_debug_msg(__func__, "Cleaning done: %.2f.\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	free(kmer_list);
	free(counter_list);
	free(sorted_counters);
	destroy_ht(ht);
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
