#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"

#ifdef HAVE_PTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

gap_opt_t *gap_init_opt() {
	gap_opt_t *o;
	o = (gap_opt_t*) calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	 rate. Voilating this requirement will break pairing! */
	o->s_mm = 3;
	o->s_gapo = 11;
	o->s_gape = 4;
	o->max_diff = -1;
	o->max_gapo = 1;
	o->max_gape = 6;
	o->indel_end_skip = 5;
	o->max_del_occ = 10;
	o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32;
	o->max_seed_diff = 2;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres) {
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres)
			return k;
	}
	return 2;
}

// width must be filled as zero
static int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str,
		bwt_width_t *width) {
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0;
	l = rbwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
			k = rbwt->L2[c] + ok + 1;
			l = rbwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = rbwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

void bwa_cal_sa_reg_gap(int tid, bwt_t * const bwt[2], int n_seqs,
		bwa_seq_t *seqs, const gap_opt_t *opt) {
	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len)
			max_len = seqs[i].len;
	if (opt->fnr > 0.0)
		local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo)
		local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo,
			local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*) calloc(opt->seed_len + 1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*) calloc(opt->seed_len + 1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
#ifdef HAVE_PTHREAD
		if (opt->n_threads > 1) {
			pthread_mutex_lock(&g_seq_lock);
			if (p->tid < 0) { // unassigned
				int j;
				for (j = i; j < n_seqs && j < i + THREAD_BLOCK_SIZE; ++j)
					seqs[j].tid = tid;
			} else if (p->tid != tid) {
				pthread_mutex_unlock(&g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&g_seq_lock);
		}
#endif
		seq[0] = p->seq;
		seq[1] = p->rseq;
		if (max_l < p->len) {
			max_l = p->len;
			w[0] = (bwt_width_t*) realloc(w[0], (max_l + 1)
					* sizeof(bwt_width_t));
			w[1] = (bwt_width_t*) realloc(w[1], (max_l + 1)
					* sizeof(bwt_width_t));
			memset(w[0], 0, (max_l + 1) * sizeof(bwt_width_t));
			memset(w[1], 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], p->len, seq[0], w[0]);
		bwt_cal_width(bwt[1], p->len, seq[1], w[1]);
		if (opt->fnr > 0.0)
			local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len ? opt->seed_len
				: 0x7fffffff;
		if (p->len > opt->seed_len) {
			bwt_cal_width(bwt[0], opt->seed_len, seq[0] + (p->len
					- opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len, seq[1] + (p->len
					- opt->seed_len), seed_w[1]);
		}
		free(p->name);
		free(p->seq);
		free(p->rseq);
		p->name = 0;
		p->seq = p->rseq = 0;
	}
	free(seed_w[0]);
	free(seed_w[1]);
	free(w[0]);
	free(w[1]);
	gap_destroy_stack(stack);
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bwt_t *bwt[2];
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data) {
	thread_aux_t *d = (thread_aux_t*) data;
	bwa_cal_sa_reg_gap(d->tid, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}
#endif

bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa) {
	bwa_seqio_t *ks;
	if (mode & BWA_MODE_BAM) { // open BAM
		int which = 0;
		if (mode & BWA_MODE_BAM_SE)
			which |= 4;
		if (mode & BWA_MODE_BAM_READ1)
			which |= 1;
		if (mode & BWA_MODE_BAM_READ2)
			which |= 2;
		if (which == 0)
			which = 7; // then read all reads
		ks = bwa_bam_open(fn_fa, which);
	} else
		ks = bwa_seq_open(fn_fa);
	return ks;
}

bwa_alg *pe_aln_core(bwa_seq_t *seqs, const gap_opt_t *opt,
		const bwt_t *bwt_0, const bwt_t *bwt_1) {
	clock_t t;
	bwt_t *bwt[2];
	bwa_alg *align = (bwa_alg*) malloc(sizeof(bwa_alg));

	bwt[0] = bwt_0;
	bwt[1] = bwt_1;
	t = clock();
	int n_seqs = 1;

#ifdef HAVE_PTHREAD
	if (opt->n_threads <= 1) { // no multi-threading at all
		bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
	} else {
		pthread_t *tid;
		pthread_attr_t attr;
		thread_aux_t *data;
		int j;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		data = (thread_aux_t*) calloc(opt->n_threads, sizeof(thread_aux_t));
		tid = (pthread_t*) calloc(opt->n_threads, sizeof(pthread_t));
		for (j = 0; j < opt->n_threads; ++j) {
			data[j].tid = j;
			data[j].bwt[0] = bwt[0];
			data[j].bwt[1] = bwt[1];
			data[j].n_seqs = n_seqs;
			data[j].seqs = seqs;
			data[j].opt = opt;
			pthread_create(&tid[j], &attr, worker, data + j);
		}
		for (j = 0; j < opt->n_threads; ++j)
			pthread_join(tid[j], 0);
		free(data);
		free(tid);
	}
#else
	bwa_cal_sa_reg_gap(0, bwt, n_seqs, seqs, opt);
#endif

	t = clock();
	bwa_free_read_seq(n_seqs, seqs);
	return align;
}

/* rgoya: Temporary clone of aln_path2cigar to accomodate for bwa_cigar_t,
 __cigar_op and __cigar_len while keeping stdaln stand alone */
bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar) {
	uint32_t *cigar32;
	bwa_cigar_t *cigar;
	int i;
	cigar32 = aln_path2cigar32((path_t*) path, path_len, n_cigar);
	cigar = (bwa_cigar_t*) cigar32;
	for (i = 0; i < *n_cigar; ++i)
		cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}

