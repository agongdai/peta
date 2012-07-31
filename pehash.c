/*
 * pehash.c
 *
 *  Created on: 11-Dec-2011
 *      Author: carl
 */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <glib.h>
#include <inttypes.h>
#include "pehash.h"
#include "bwase.h"
#include "peseq.h"

hash_opt *init_hash_opt() {
	hash_opt *o = (hash_opt*) calloc(1, sizeof(hash_opt));
	o->mode = (BWA_MODE_GAPE | BWA_MODE_COMPREAD);
	o->k = 11;
	o->read_len = 0;
	o->n_k_mers = 0;
	o->n_pos = 0;
	o->interleaving = 0;
	o->n_hash_block = 2;
	o->block_size = o->k / 2;
	return o;
}

int hash_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta hash [options] <rnaseq.fa> \n");
	return 1;
}

hash_key get_hash_key(const ubyte_t *seq, const int start,
		const int interleaving, const int len) {
	hash_key key = 0ULL;
	int i = 0;
	for (i = 0; i < len * interleaving; i += interleaving) {
		if (seq[i + start] < 0 || seq[i + start] > 3) {
			return 0;
		}
		key <<= 2;
		// printf("%d \n", seq[i + start]);
		key = key | seq[i + start];
		// printf("key = %" ID64 "\n", key);
	}
	return key;
}

/**
 * The elements in pos array are 64bits, where upper 48 bits representing the rna-seq id,
 * 	and lower 16 bits representing the starting position of hashing sub-sequence.
 */
hash_value get_hash_value(const index64 seq_id, const int pos_start) {
	hash_value value = seq_id;
	value <<= N_POS_BITS;
	value = value | pos_start;
	// printf("[get_hash_value] (%" ID64 ", %d) => %" ID64 "\n", seq_id, pos_start, value);
	return value;
}

void read_hash_value(index64 *seq_id, int *pos_start, hash_value value) {
	*seq_id = value >> N_POS_BITS;
	*seq_id = *seq_id & HASH_VALUE_HIGHER;
	*pos_start = value & HASH_VALUE_LOWER;
}

void pe_hash_core(const char *fa_fn, hash_opt *opt) {
	bwa_seq_t *s, *part_seqs;
	bwa_seqio_t *ks;
	index64 key, n_k_mers = 0;
	hash_value value;
	int hash_len = 0, tmp = 0, tmp_2 = 0, hash_start = 0, block_no = 0;
	index64 i = 0, i_acc = 0, n_pos = 0, pos_index = 0;
	uint64_t n_part_seqs = 0, n_seqs = 0;
	clock_t t = clock();
	FILE *hash_fp;
	hash_key *k_mers, *k_mers_occ_acc;
	hash_value *pos;
	char *hash_fn = malloc(BUFSIZE);
	int *n_occ;

	// All k-mer combinations
	n_k_mers = (1 << (opt->k * 2)) + 1;
	fprintf(stderr,
			"[pe_hash_core] K-mer size: %d; k-mer hashtable size: %" ID64 "\n",
			opt->k, n_k_mers);
	k_mers = (hash_value*) calloc(n_k_mers, sizeof(hash_value));
	k_mers_occ_acc = (hash_value*) calloc(n_k_mers, sizeof(hash_value));
	n_occ = (int*) calloc(n_k_mers, sizeof(int)); // Store how many occs through the iteration
	hash_len = (opt->read_len / opt->k) * opt->k;

	fprintf(stderr, "[pe_hash_core] Loading sequences from %s... \n", fa_fn);
	ks = bwa_open_reads(opt->mode, fa_fn);
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, opt->mode,
			0)) != 0) {
		pe_reverse_seqs(part_seqs, n_part_seqs);
		n_seqs += n_part_seqs;
		fprintf(
				stderr,
				"[pe_hash_core] %"ID64" sequences in library %s loaded: %.2f sec... \n",
				n_seqs, fa_fn, (float) (clock() - t)
				/ CLOCKS_PER_SEC);
		//		p_query("", part_seqs);
		// Round 1: count occurrences of all k-mers
		fprintf(stderr,
				"[pe_hash_core] Round 1: Counting occurrences of k-mers... \n");
		for (i = 0; i < n_part_seqs; i++) {
			s = &part_seqs[i];
			if (s->len != opt->read_len) {
				err_fatal(
						__func__,
						"Sequence length of %s is not as specified: %d vs %d! \n",
						s->name, s->len, opt->read_len);
			}
			//p_query(__func__, s);
			while (block_no < opt->n_hash_block && hash_start <= opt->read_len
					- opt->k * (opt->interleaving)) {
				hash_start = block_no * opt->block_size;
				key = get_hash_key(s->seq, hash_start, opt->interleaving,
						opt->k);
				//printf("Key: %" ID64 "\n", key);
				k_mers_occ_acc[key]++;
				n_pos++;

				key = get_hash_key(s->seq, hash_start + opt->interleaving - 1,
						opt->interleaving, opt->k);
				k_mers_occ_acc[key]++;
				n_pos++;
				block_no++;
			}
			hash_start = 0;
			block_no = 0;
		}
		bwa_free_read_seq(n_part_seqs, part_seqs);
	}
	bwa_seq_close(ks);
	n_seqs = 0;

	// Example:
	//	k-mers: 1, 2, 3, 4;
	//	occurrences for each k-mer: 2, 4, 6, 8;
	//	accumulate all occs: 0, 2, 6, 12 (end with the last index)
	tmp = 0;
	for (i = 1; i < n_k_mers - 1; i++) {
		tmp_2 = k_mers_occ_acc[i];
		k_mers_occ_acc[i] = k_mers_occ_acc[i - 1] + tmp;
		tmp = tmp_2;
	}
	k_mers_occ_acc[0] = 0;
	k_mers_occ_acc[n_k_mers - 1] = n_pos - 1;

	// Round 2: save the hashes.
	fprintf(stderr,
			"[pe_hash_core] Pos array size: %" ID64 ", start hashing...\n",
			n_pos);
	pos = (hash_value*) calloc(n_pos, sizeof(hash_value));

	hash_start = 0;
	block_no = 0;
	fprintf(stderr, "[pe_hash_core] Round 2: Loading sequences from %s... \n",
			fa_fn);
	ks = bwa_open_reads(opt->mode, fa_fn);
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, opt->mode,
			0)) != 0) {
		pe_reverse_seqs(part_seqs, n_part_seqs);
		n_seqs += n_part_seqs;
		fprintf(
				stderr,
				"[pe_hash_core] %"ID64" sequences in library %s loaded: %.2f sec... \n",
				n_seqs, fa_fn, (float) (clock() - t)
				/ CLOCKS_PER_SEC);

		for (i = 0; i < n_part_seqs; i++) {
			s = &part_seqs[i];
			while (block_no < opt->n_hash_block && hash_start <= opt->read_len
					- opt->k * (opt->interleaving)) {
				hash_start = block_no * opt->block_size;
				// Hash two keys for the starting region of a read, interleaving by 1 by default.
				key = get_hash_key(s->seq, hash_start, opt->interleaving,
						opt->k);
				value = get_hash_value(i_acc, hash_start);
				pos_index = k_mers_occ_acc[key] + n_occ[key];
				pos[pos_index] = value;
				n_occ[key]++;

				key = get_hash_key(s->seq, hash_start + opt->interleaving - 1,
						opt->interleaving, opt->k);
				value = get_hash_value(i_acc, hash_start + opt->interleaving - 1);
				pos_index = k_mers_occ_acc[key] + n_occ[key];
				pos[pos_index] = value;
				n_occ[key]++;

				block_no++;
			}
			hash_start = 0;
			block_no = 0;
			i_acc++;
		}
		bwa_free_read_seq(n_part_seqs, part_seqs);
	}
	bwa_seq_close(ks);

	sprintf(hash_fn, "%s.hash", fa_fn);
	hash_fp = xopen(hash_fn, "w");
	opt->n_k_mers = n_k_mers;
	opt->n_pos = n_pos;
	fwrite(opt, sizeof(hash_opt), 1, hash_fp);
	fprintf(stderr, "[pe_hash_core] Saving k-mers... \n");
	fwrite(k_mers_occ_acc, sizeof(hash_key), n_k_mers, hash_fp);
	fprintf(stderr, "[pe_hash_core] Saving occurrences... \n");
	fwrite(pos, sizeof(hash_value), n_pos, hash_fp);

	fclose(hash_fp);
	free(k_mers_occ_acc);
	free(n_occ);
}

hash_table *pe_load_hash(const char *fa_fn) {
	hash_table *h;
	FILE *fp;
	hash_opt *opt;
	bwa_seqio_t *ks;
	uint32_t n_seqs = 0, n_seqs_full = 0;
	int n_part_seqs = 0;
	bwa_seq_t *seqs = NULL, *part_seqs = NULL;
	char *hash_fn = malloc(FNLEN);
	clock_t t = clock();

	fprintf(stderr, "[pe_load_hash] Loading hash table of %s... \n", fa_fn);
	h = (hash_table*) malloc(sizeof(hash_table));
	// Load original rna-seqs.
	ks = bwa_open_reads(BWA_MODE, fa_fn);
	n_seqs_full = N_DF_MAX_SEQS;
	seqs = (bwa_seq_t*) calloc(N_DF_MAX_SEQS, sizeof(bwa_seq_t));
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, BWA_MODE,
			0)) != 0) {
		pe_reverse_seqs(part_seqs, n_part_seqs);
		fprintf(
				stderr,
				"[pe_load_hash] %d sequences in library %s loaded: %.2f sec... \n",
				n_seqs + n_part_seqs, fa_fn,
				(float) (clock() - t) / CLOCKS_PER_SEC);
		if ((n_seqs + n_part_seqs) > n_seqs_full) {
			n_seqs_full += n_part_seqs + 2;
			kroundup32(n_seqs_full);
			seqs = (bwa_seq_t*) realloc(seqs, sizeof(bwa_seq_t) * n_seqs_full);
		}
		memmove(&seqs[n_seqs], part_seqs, sizeof(bwa_seq_t) * n_part_seqs);
		free(part_seqs);
		n_seqs += n_part_seqs;
	}
	h->seqs = seqs;
	h->n_seqs = n_seqs;
	fprintf(stderr, "[pe_load_hash] Original sequences loaded: %d\n", n_seqs);

	// Load the hash table itself
	sprintf(hash_fn, "%s.hash", fa_fn);
	opt = (hash_opt*) malloc(sizeof(hash_opt));
	fp = xopen(hash_fn, "rb");
	if (!fread(opt, sizeof(hash_opt), 1, fp)) {
		err_fatal(__func__, "Unable to read from the hash file %s! \n", hash_fn);
	}fprintf(
			stderr,
			"[pe_load_hash] Hashing options: k=%d, read_len=%d, n_k_mers=%" ID64 ", n_pos=%" ID64 "...\n",
			opt->k, opt->read_len, opt->n_k_mers, opt->n_pos);
	h->k_mers_occ_acc = (hash_key*) calloc(opt->n_k_mers, sizeof(hash_key));
	h->pos = (hash_value*) calloc(opt->n_pos, sizeof(hash_value));
	h->o = opt;
	fread(h->k_mers_occ_acc, sizeof(hash_key), opt->n_k_mers, fp);
	fread(h->pos, sizeof(hash_value), opt->n_pos, fp);
	fclose(fp);
	fprintf(
			stderr,
			"[pe_load_hash] Hash table loaded, k-mer records: %" ID64 ", positions: %" ID64 "\n",
			opt->n_k_mers, opt->n_pos);
	bwa_seq_close(ks);
	free(hash_fn);
	return h;
}

void pe_hash_test(const char *fa_fn, hash_opt *opt) {
	hash_table *ht = NULL;
	int i = 0;
	bwa_seq_t *s;
	hash_opt *o;

	show_debug_msg(__func__, "=======================================\n");
	ht = pe_load_hash(fa_fn);
	o = ht->o;

	for (i = 0; i < o->n_k_mers; i++) {
		if (i % (200000) == 0) {
			printf("%d: %d\n", i, ht->k_mers_occ_acc[i]);
		}
	}
	for (i = 0; i < ht->n_seqs; i++) {
		if (i % (N_CHUNK_SEQS / 10) == 0) {
			s = &ht->seqs[i];
			p_query(__func__, s);
		}
	}
	for (i = 0; i < o->n_pos; i++) {
		if (i % (N_CHUNK_SEQS) == 0) {
			printf("%d: %d\n", i, ht->pos[i]);
		}
	}
}

int pe_hash(int argc, char *argv[]) {
	int c;
	clock_t t;
	t = clock();
	hash_opt *opt = init_hash_opt();

	while ((c = getopt(argc, argv, "k:l:i:b:s:")) >= 0) {
		switch (c) {
		case 'k':
			opt->k = atoi(optarg);
			break;
		case 'l':
			opt->read_len = atoi(optarg);
			break;
		case 'i':
			opt->interleaving = atoi(optarg);
			break;
		case 'b':
			opt->n_hash_block = atoi(optarg);
			break;
		case 's':
			opt->block_size = atoi(optarg);
			break;
		default:
			return 1;
		}
	}

	if (optind + 1 > argc) {
		return hash_usage();
	}

	pe_hash_core(argv[optind], opt);
//	pe_hash_test(argv[optind], opt);
//	pe_hash_test("read/SRR097897.fa.correct", opt);
	fprintf(stderr, "[pe_hash] Hashing done: %.2f sec\n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return 0;
}
