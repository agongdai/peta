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
#include "rnaseq.h"

GMutex *counter_mutex;

int clean_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta clean [options] \n");
	return 1;
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
	o->lib_name = NULL;
	o->mode = (BWA_MODE_GAPE | BWA_MODE_COMPREAD);
	o->n_threads = 1;
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

/**
 * Check whether some k-mer in the read is extremely rare (fewer than 2).
 */
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

/**
 * A read 80bp, k-mer length 10.
 * For each base in the read, count the k-mer frequency in the middle [10, 60]
 * If a k-mer is [20, 30], count the base k-mer frequencies of each base
 */
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

/**
 * Calculate the mean and standard deviation of the k-mer frequency
 */
void set_k_freq(bwa_seq_t *read, counter *k_count, uint16_t *kmer_list,
		const int k) {
	hash_key key = 0;
	int i = 0;
	double counted_len = read->len - k;
	double *base_counter = (double*) calloc(read->len - k + 1, sizeof(double));
	for (i = 0; i <= read->len - k; i++) {
		key = get_key(read, i, i + k);
		k_count->k_freq += kmer_list[key];
	}
	for (i = 0; i <= read->len - k; i++) {
		key = get_key(read, i, i + k);
		base_counter[i] = kmer_list[key];
	}
	k_count->k_mean = mean(base_counter, counted_len);
	k_count->k_sd = std_dev(base_counter, counted_len);
	free(base_counter);
}

void *set_k_freq_t(void *data) {
	clean_thread_aux *d = (clean_thread_aux*) data;
	bwa_seq_t *seqs = d->seqs, *s = NULL;
	int i = 0;
	counter *k_count = NULL, *counter_list = d->counter_list;

	for (i = d->start; i < d->end; i++) {
		s = &seqs[i];
		if (s->status != FRESH && s->status != TRIED)
			continue;
		k_count = &counter_list[i];
		if (k_count->checked) {
			continue;
		}
		set_k_freq(s, k_count, d->kmer_list, d->opt->kmer);
	}
	return NULL;
}

void set_k_freq_threads(bwa_seq_t *seqs, const int n_seqs,
		counter *counter_list, uint16_t *kmer_list, clean_opt *opt) {
	int i = 0, n_per_threads = 0;
	GThread *threads[opt->n_threads];
	clean_thread_aux *data;

	data = (clean_thread_aux*) calloc(opt->n_threads, sizeof(clean_thread_aux));
	n_per_threads = n_seqs / opt->n_threads;
	for (i = 0; i < opt->n_threads; ++i) {
		data[i].seqs = seqs;
		data[i].n_seqs = n_seqs;
		data[i].counter_list = counter_list;
		data[i].kmer_list = kmer_list;
		data[i].opt = opt;
		data[i].start = n_per_threads * i;
		data[i].end = n_per_threads * (i + 1);
		if (i == opt->n_threads - 1)
			data[i].end = n_seqs;
		threads[i] = g_thread_create((GThreadFunc) set_k_freq_t, data + i,
				TRUE, NULL);
	}
	for (i = 0; i < opt->n_threads; ++i) {
		g_thread_join(threads[i]);
	}
	free(data);
}

void *rm_rep_low(void *data) {
	clean_thread_aux *d = (clean_thread_aux*) data;
	bwa_seq_t *seqs = d->seqs, *s = NULL;
	int i = 0;
	counter *k_count, *counter_list = d->counter_list;
	for (i = d->start; i < d->end; i++) {
		s = &seqs[i];
		if (s->status != FRESH && s->status != TRIED)
			continue;
		k_count = &counter_list[i];
		if (k_count->checked) {
			continue;
		}
		if (is_biased_q(s) || has_rep_pattern(s) || is_repetitive_q(s)) {
			k_count->checked = 1;
		} else {
			if (d->rm_low_kmer && has_low_kmer(s, d->kmer_list, d->opt)) {
				k_count->checked = 2;
			}
		}
	}
	return NULL;
}

void rm_rep_low_threads(bwa_seq_t *seqs, const int n_seqs, clean_opt *opt,
		uint16_t *kmer_list, counter *counter_list, const int rm_low_kmer) {
	int i = 0, n_per_threads = 0;
	GThread *threads[opt->n_threads];
	clean_thread_aux *data;

	data = (clean_thread_aux*) calloc(opt->n_threads, sizeof(clean_thread_aux));
	n_per_threads = n_seqs / opt->n_threads;
	for (i = 0; i < opt->n_threads; ++i) {
		data[i].seqs = seqs;
		data[i].n_seqs = n_seqs;
		data[i].rm_low_kmer = rm_low_kmer;
		data[i].counter_list = counter_list;
		data[i].kmer_list = kmer_list;
		data[i].opt = opt;
		data[i].start = n_per_threads * i;
		data[i].end = n_per_threads * (i + 1);
		if (i == opt->n_threads - 1)
			data[i].end = n_seqs;
		threads[i] = g_thread_create((GThreadFunc) rm_rep_low, data + i, TRUE,
				NULL);
	}
	for (i = 0; i < opt->n_threads; ++i) {
		g_thread_join(threads[i]);
	}
	free(data);
}

void *rm_dup(void *data) {
	int i = 0;
	iterate_thread_aux *d = (iterate_thread_aux*) data;
	bwa_seq_t *seqs = d->seqs, *s = NULL, *s_unique = NULL;
	counter *sorted_counters = d->sorted_counters, *k_count = NULL,
			*counter_pre = sorted_counters;
	s_unique = &seqs[d->start];
	for (i = d->start; i < d->end; i++) {
		k_count = &sorted_counters[i];
		s = &seqs[k_count->read_id];
		if (s->status != FRESH && s->status != TRIED)
			continue;
		if (k_count->k_freq == counter_pre->k_freq && same_q(s, s_unique)) {
			k_count->checked = 3;
		} else {
			s_unique = s;
		}
		counter_pre = k_count;
	}
}

void rm_dup_threads(bwa_seq_t *seqs, const int n_seqs, clean_opt *opt,
		counter *sorted_counters) {
	int i = 0, n_per_threads = 0;
	GThread *threads[opt->n_threads];
	iterate_thread_aux *data;

	data = (iterate_thread_aux*) calloc(opt->n_threads,
			sizeof(iterate_thread_aux));
	n_per_threads = n_seqs / opt->n_threads;
	for (i = 0; i < opt->n_threads; ++i) {
		data[i].seqs = seqs;
		data[i].n_seqs = n_seqs;
		data[i].sorted_counters = sorted_counters;
		data[i].start = n_per_threads * i;
		data[i].end = n_per_threads * (i + 1);
		if (i == opt->n_threads - 1)
			data[i].end = n_seqs;
		threads[i]
				= g_thread_create((GThreadFunc) rm_dup, data + i, TRUE, NULL);
	}
	for (i = 0; i < opt->n_threads; ++i) {
		g_thread_join(threads[i]);
	}
	free(data);
}

void *iterate_seqs(void *data) {
	iterate_thread_aux *d = (iterate_thread_aux*) data;
	int i = 0, *n_solid = d->n_solid;
	counter *sorted_counters = d->sorted_counters, *k_count = NULL;
	bwa_seq_t *seqs = d->seqs, *s = NULL;
	for (i = d->start; i < d->end; i++) {
		k_count = &sorted_counters[i];
		if (k_count->checked)
			continue;
		s = &seqs[k_count->read_id];
		if (s->status != FRESH && s->status != TRIED)
			continue;
		if (pick_within_range(s, k_count, d->kmer_list, d->opt, UNEVEN_THRE
				* d->n_iterate)) {
			k_count->checked = 4;
			g_mutex_lock(counter_mutex);
			*n_solid += 1;
			g_ptr_array_add(d->solid_reads, s);
			g_mutex_unlock(counter_mutex);
		}
		if (*n_solid >= d->n_needed)
			break;
	}
}

void iterate_seqs_threads(bwa_seq_t *seqs, const int n_seqs, clean_opt *opt,
		counter *sorted_counters, uint16_t *kmer_list, int *n_solid,
		int n_needed, GPtrArray *solid_reads, const int n_iterate) {
	int i = 0, n_per_threads = 0;
	GThread *threads[opt->n_threads];
	iterate_thread_aux *data;

	data = (iterate_thread_aux*) calloc(opt->n_threads,
			sizeof(iterate_thread_aux));
	n_per_threads = n_seqs / opt->n_threads;
	for (i = 0; i < opt->n_threads; ++i) {
		data[i].seqs = seqs;
		data[i].n_seqs = n_seqs;
		data[i].n_solid = n_solid;
		data[i].n_needed = n_needed;
		data[i].sorted_counters = sorted_counters;
		data[i].kmer_list = kmer_list;
		data[i].opt = opt;
		data[i].start = n_per_threads * i;
		data[i].end = n_per_threads * (i + 1);
		data[i].n_iterate = n_iterate;
		data[i].solid_reads = solid_reads;
		if (i == opt->n_threads - 1)
			data[i].end = n_seqs;
		threads[i] = g_thread_create((GThreadFunc) iterate_seqs, data + i,
				TRUE, NULL);
	}
	for (i = 0; i < opt->n_threads; ++i) {
		g_thread_join(threads[i]);
	}
	free(data);
}

GPtrArray *calc_solid_reads(bwa_seq_t *seqs, const int n_seqs, clean_opt *opt,
		const int n_needed, const int by_coverage, const int rm_low_kmer) {
	int i = 0, j = 0;
	int n_solid = 0, n_has_n = 0, try_times = 0;
	uint32_t n_kmers = 0;
	uint16_t *kmer_list;
	counter *k_count = NULL, *counter_list = NULL, *sorted_counters = NULL;
	GPtrArray *solid_reads = NULL;
	bwa_seq_t *s = NULL;
	struct timespec start, finish;
	// Each counter corresponds to a read.
	// It contains the kmer frequencies.
	counter_list = (counter*) calloc(n_seqs, sizeof(counter));

	if (!g_thread_supported())
		g_thread_init(NULL);
	if (!counter_mutex)
		counter_mutex = g_mutex_new();
	clock_gettime(CLOCK_MONOTONIC, &start);

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
	clock_gettime(CLOCK_MONOTONIC, &finish);
	show_debug_msg(__func__, "Calculating k-mer frequency: %.2f sec...\n",
			(float) (finish.tv_sec - start.tv_sec));
	for (i = 0; i < n_seqs; i++) {
		s = &seqs[i];
		if (s->status != FRESH && s->status != TRIED) {
			continue;
		}
		k_count = &counter_list[i];
		k_count->read_id = atoi(s->name);
		k_count->checked = 0;
		if (has_n(s)) {
			k_count->checked = -1;
			n_has_n++;
			continue;
		}
		set_kmer_index(s, opt->kmer, kmer_list);
	}

	show_debug_msg(__func__, "At most %d reads are being targeted. \n",
			n_needed);
	clock_gettime(CLOCK_MONOTONIC, &finish);
	show_debug_msg(__func__, "%d reads with 'N' removed: %.2f sec...\n",
			n_has_n, (float) (finish.tv_sec - start.tv_sec));

	// Calculate the mean and standard deviation of k_freq
	show_debug_msg(__func__,
			"Counting the k-mer frequencies of reads: %.2f sec...\n",
			(float) (finish.tv_sec - start.tv_sec));
	set_k_freq_threads(seqs, n_seqs, counter_list, kmer_list, opt);

	clock_gettime(CLOCK_MONOTONIC, &finish);
	show_debug_msg(__func__,
			"Removing repetitive and low k-mer frequency reads: %.2f sec...\n",
			(float) (finish.tv_sec - start.tv_sec));
	// Remove repetitive reads and those reads having low frequency kmers.
	rm_rep_low_threads(seqs, n_seqs, opt, kmer_list, counter_list, rm_low_kmer);

	// Sort the counters
	clock_gettime(CLOCK_MONOTONIC, &finish);
	show_debug_msg(__func__, "Sorting reads by k_freq: %.2f sec...\n",
			(float) (finish.tv_sec - start.tv_sec));
	qsort(counter_list, n_seqs, sizeof(counter), cmp_kmer);
	sorted_counters = counter_list;

	clock_gettime(CLOCK_MONOTONIC, &finish);
	show_debug_msg(__func__, "Removing duplicates: %.2f sec...\n",
			(float) (finish.tv_sec - start.tv_sec));
	rm_dup_threads(seqs, n_seqs, opt, sorted_counters);
	clock_gettime(CLOCK_MONOTONIC, &finish);

	show_debug_msg(__func__, "Getting solid reads: %.2f sec...\n",
			(float) (finish.tv_sec - start.tv_sec));
	solid_reads = g_ptr_array_sized_new(16384);
	if (by_coverage) {
		for (i = 0; i < n_seqs; i++) {
			k_count = &sorted_counters[i];
			if (k_count->checked)
				continue;
			s = &seqs[k_count->read_id];
			if (s->status != FRESH && s->status != TRIED)
				continue;
			n_solid += 1;
			g_ptr_array_add(solid_reads, s);
			if (n_solid >= n_needed)
				break;
		}
	} else {
		j = 0;
		try_times = MAX_TIME;
		while (++j <= try_times && n_needed >= n_solid) {
			// For low sd range, there are only few reads are solid
			// Here is to avoid unnecessary loops on the reads.
			clock_gettime(CLOCK_MONOTONIC, &finish);
			show_debug_msg(__func__,
					"Round %d out of max %d. %d solid reads: %.2f sec...\n", j,
					try_times, n_solid, (float) (finish.tv_sec - start.tv_sec));
			iterate_seqs_threads(seqs, n_seqs, opt, sorted_counters, kmer_list,
					&n_solid, n_needed, solid_reads, j);
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &finish);
	show_debug_msg(__func__, "%d solid reads remained.\n", n_solid);
	show_debug_msg(__func__, "Cleaning done: %.2f.\n", (float) (finish.tv_sec
			- start.tv_sec));
	free(kmer_list);
	free(counter_list);
	return solid_reads;
}

void pe_clean_core(char *fa_fn, clean_opt *opt) {
	bwa_seq_t *seqs, *s = NULL;
	int n_seqs = 0, i = 0;
	char *item = (char*) malloc(BUFSIZE), *solid = malloc(BUFSIZE);
	FILE *solid_file;
	clock_t t = clock();
	GPtrArray *solid_reads = NULL;

	show_debug_msg(__func__, "Loading library %s...\n", fa_fn);
	seqs = load_reads(fa_fn, &n_seqs);

	show_debug_msg(__func__, "Saving k-mer frequencies: %.2f sec...\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	sprintf(solid, "%s.solid", opt->lib_name);
	solid_file = xopen(solid, "w");
	solid_reads = calc_solid_reads(seqs, n_seqs, opt, n_seqs * opt->stop_thre,
			0, 1);
	for (i = 0; i < solid_reads->len; i++) {
		s = g_ptr_array_index(solid_reads, i);
		sprintf(item, "%s\n", s->name);
		fputs(item, solid_file);
	}

	free(item);
	free(solid);
	g_ptr_array_free(solid_reads, TRUE);
	bwa_free_read_seq(n_seqs, seqs);
	fclose(solid_file);
}

int clean_reads(int argc, char *argv[]) {
	int c;
	struct timespec start, finish;
	clean_opt *opt = init_clean_opt();

	clock_gettime(CLOCK_MONOTONIC, &start);
	while ((c = getopt(argc, argv, "k:l:s:t:")) >= 0) {
		switch (c) {
		case 'k':
			opt->kmer = atoi(optarg);
			break;
		case 'l':
			opt->lib_name = optarg;
			break;
		case 's':
			opt->stop_thre = atof(optarg);
			break;
		case 't':
			opt->n_threads = atoi(optarg);
			break;
		default:
			return 1;
		}
	}
	if (!g_thread_supported())
		g_thread_init(NULL);

	if (optind + 1 > argc) {
		return clean_usage();
	}
	//test_threads();
	pe_clean_core(argv[optind], opt);
	free(opt);
	clock_gettime(CLOCK_MONOTONIC, &finish);
	fprintf(stderr, "[clean_reads] Cleaning done: %.2f sec\n",
			(float) (finish.tv_sec - start.tv_sec));
	return 0;
}
