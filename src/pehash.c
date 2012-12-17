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
#include "rnaseq.h"
#include "edgelist.h"

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

void destroy_ht(hash_table *ht) {
	if (ht) {
		free(ht->o);
		free(ht->k_mers_occ_acc);
		free(ht->pos);
		bwa_free_read_seq(ht->n_seqs, ht->seqs);
		free(ht);
	}
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

void g_array_remove_fast(GArray *arr, hash_value value) {
	int v = 0, i = 0;
	for (i = 0; i < arr->len; i++) {
		v = g_array_index(arr, hash_value, i);
		if (v == value) {
			if (g_array_remove_index_fast(arr, i))
				i--;
		}
	}
}

void add_read_to_ht(reads_ht *ht, bwa_seq_t *read) {
	int i = 0;
	hash_key key = 0;
	hash_value value = 0;
	GArray *occs = NULL;

	//p_query(__func__, read);
	for (i = 0; i < read->len - ht->k; i++) {
		key = get_hash_key(read->seq, i, 1, ht->k);
		value = get_hash_value(atoi(read->name), i);
		//show_debug_msg(__func__, "ADDED ENTRY: %" ID64 "=>%" ID64 " \n", key,
		//		value);
		occs = g_ptr_array_index(ht->pos, key);
		g_array_append_val(occs, value);
	}
	ht->n_reads++;
}

void rm_read_from_ht(reads_ht *ht, bwa_seq_t *read) {
	int i = 0;
	hash_key key = 0;
	hash_value value = 0;
	GArray *occs = NULL;
	if (!ht)
		return;
	for (i = 0; i < read->len - ht->k; i++) {
		key = get_hash_key(read->seq, i, 1, ht->k);
		value = get_hash_value(atoi(read->name), i);
		occs = g_ptr_array_index(ht->pos, key);
		g_array_remove_fast(occs, value);
	}
	ht->n_reads--;
}

GPtrArray *mutate_one_base(ubyte_t *seq, const int start_index, const int len) {
	GPtrArray *mutated_all = NULL;
	ubyte_t *mutated = NULL, *target_sub_seq = NULL;
	int j = 0, i = 0;
	mutated_all = g_ptr_array_sized_new(len * 3 + 1);
	target_sub_seq = (ubyte_t*) calloc(len + 1, sizeof(ubyte_t));
	memcpy(target_sub_seq, &seq[start_index], len);
	target_sub_seq[len] = '\0';
	//p_seq(__func__, target_sub_seq, len);
	for (i = 0; i < len; i++) {
		// i: the position to mutate;
		// j: four possible nucleotides: 'a'=>0, 'c'=>1, ...
		for (j = 0; j < 4; j++) {
			if (target_sub_seq[i] != j) {
				mutated = (ubyte_t*) calloc(len + 1, sizeof(ubyte_t));
				memcpy(mutated, target_sub_seq, len);
				mutated[len] = '\0';
				mutated[i] = j;
				g_ptr_array_add(mutated_all, mutated);
			}
		}
	}
	g_ptr_array_add(mutated_all, target_sub_seq);
	return mutated_all;
}

GPtrArray *find_ol_reads(reads_ht *ht, bwa_seq_t *template, bwa_seq_t *seqs,
		GPtrArray *hit_reads, const int ori) {
	hash_key key = 0;
	GArray *occs = NULL;
	int i = 0, j = 0, locus = 0;
	index64 seq_id = 0;
	hash_value value = 0;
	bwa_seq_t *r = NULL;
	GPtrArray *mutated = NULL;
	ubyte_t *tmp = NULL, *m = NULL;
	//show_debug_msg(__func__, "Getting mutated template... \n");
	if (ori) {
		mutated = mutate_one_base(template->seq, 0, ht->k);
		key = get_hash_key(template->seq, 0, 1, ht->k);
	} else {
		mutated = mutate_one_base(template->seq, template->len - ht->k, ht->k);
		key = get_hash_key(template->seq, template->len - ht->k, 1, ht->k);
	}
	for (j = 0; j < ht->k * 3 + 1; j++) {
		if (j == 0) {
			//p_seq(__func__, template->seq, template->len);
			//show_debug_msg(__func__, "KEY: %" ID64 " \n", key);
		} else {
			m = g_ptr_array_index(mutated, j - 1);
			//p_seq(__func__, m, ht->k);
			key = get_hash_key(m, 0, 1, ht->k);
			//show_debug_msg(__func__, "KEY: %" ID64 " \n", key);
		}
		occs = g_ptr_array_index(ht->pos, key);
		if (occs) {
			if (!hit_reads)
				hit_reads = g_ptr_array_sized_new(occs->len);
			for (i = 0; i < occs->len; i++) {
				value = g_array_index(occs, hash_value, i);
				//show_debug_msg(__func__, "VALUE: %" ID64 "\n", value);
				read_hash_value(&seq_id, &locus, value);
				r = &seqs[seq_id];
				//p_query("HIT", r);
				g_ptr_array_add(hit_reads, r);
			}
		}
	}
	//show_debug_msg(__func__, "Freeing mutated template... \n");
	for (i = 0; i < mutated->len; i++) {
		tmp = g_ptr_array_index(mutated, i);
		free(tmp);
	}
	g_ptr_array_free(mutated, TRUE);
	return hit_reads;
}

GPtrArray *find_reads_ol_template(reads_ht *ht, bwa_seq_t *template,
		bwa_seq_t *seqs, const int ori) {
	bwa_seq_t *rev = NULL;
	GPtrArray *hit_reads = NULL;
	if (ht->n_reads == 0) {
		hit_reads = g_ptr_array_sized_new(0);
	} else {
		//p_ctg_seq(__func__, template);
		rev = new_mem_rev_seq(template, template->len, 0);
		hit_reads = find_ol_reads(ht, template, seqs, hit_reads, ori);
		hit_reads = find_ol_reads(ht, rev, seqs, hit_reads, ori);
		bwa_free_read_seq(1, rev);
	}
	return hit_reads;
}

void destroy_reads_ht(reads_ht *ht) {
	int i = 0, n_k_mers = 0;
	GArray *occs = NULL;
	if (ht) {
		n_k_mers = (1 << (ht->k * 2)) + 1;
		for (i = 0; i < n_k_mers; i++) {
			occs = g_ptr_array_index(ht->pos, i);
			if (occs)
				g_array_free(occs, TRUE);
		}
		g_ptr_array_free(ht->pos, TRUE);
		free(ht);
	}
}

reads_ht *build_reads_ht(const int k, GPtrArray *initial_reads) {
	reads_ht *ht = NULL;
	int n_k_mers = 0, i = 0;
	GArray *occs = NULL;
	// For each k-mer, there is an arraylist storing the reads with this k-mer.
	GPtrArray *pos = NULL;
	bwa_seq_t *r = NULL;

	//show_msg(__func__, "Building hash table for mates... \n");
	n_k_mers = (1 << (k * 2)) + 1;
	ht = (reads_ht*) malloc(sizeof(reads_ht));
	pos = g_ptr_array_sized_new(n_k_mers);
	for (i = 0; i < n_k_mers; i++) {
		occs = g_array_sized_new(TRUE, TRUE, sizeof(hash_value), 0);
		g_ptr_array_add(pos, occs);
	}
	ht->pos = pos;
	ht->k = k;
	ht->n_reads = 0;
	if (initial_reads && initial_reads->len > 0) {
		for (i = 0; i < initial_reads->len; i++) {
			r = g_ptr_array_index(initial_reads, i);
			add_read_to_ht(ht, r);
		}
	}
	return ht;
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
	fprintf(stderr, "Hashing library %s...\n", fa_fn);
	fprintf(stderr,
			"[pe_hash_core] K-mer size: %d; k-mer hashtable size: %" ID64 "\n",
			opt->k, n_k_mers);
	k_mers = (hash_value*) calloc(n_k_mers, sizeof(hash_value));
	k_mers_occ_acc = (hash_value*) calloc(n_k_mers, sizeof(hash_value));
	n_occ = (int*) calloc(n_k_mers, sizeof(int)); // Store how many occs through the iteration
	hash_len = (opt->read_len / opt->k) * opt->k;

	ks = bwa_open_reads(opt->mode, fa_fn);
	// Round 1: count occurrences of all k-mers
	fprintf(stderr,
			"[pe_hash_core] Round 1/2: Counting occurrences of k-mers... \n");
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, opt->mode,
			0)) != 0) {
		pe_reverse_seqs(part_seqs, n_part_seqs);
		n_seqs += n_part_seqs;
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
		fprintf(stderr,
				"[pe_hash_core] %"ID64" sequences counted: %.2f sec... \n",
				n_seqs, (float) (clock() - t) / CLOCKS_PER_SEC);
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
	fprintf(stderr, "[pe_hash_core] Round 2/2: Store k-mer pointers %s... \n",
			fa_fn);
	ks = bwa_open_reads(opt->mode, fa_fn);
	while ((part_seqs = bwa_read_seq(ks, N_CHUNK_SEQS, &n_part_seqs, opt->mode,
			0)) != 0) {
		pe_reverse_seqs(part_seqs, n_part_seqs);
		n_seqs += n_part_seqs;
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
				value = get_hash_value(i_acc, hash_start + opt->interleaving
						- 1);
				pos_index = k_mers_occ_acc[key] + n_occ[key];
				pos[pos_index] = value;
				n_occ[key]++;

				block_no++;
			}
			hash_start = 0;
			block_no = 0;
			i_acc++;
		}
		fprintf(stderr,
				"[pe_hash_core] %"ID64" sequences hashed: %.2f sec... \n",
				n_seqs, (float) (clock() - t) / CLOCKS_PER_SEC);
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
	uint32_t n_seqs = 0;
	bwa_seq_t *seqs = NULL;
	char *hash_fn = malloc(FNLEN);
	clock_t t = clock();

	fprintf(stderr, "[pe_load_hash] Loading hash table of %s... \n", fa_fn);
	h = (hash_table*) malloc(sizeof(hash_table));
	seqs = load_reads(fa_fn, &n_seqs);
	h->seqs = seqs;
	h->n_seqs = n_seqs;
	fprintf(stderr,
			"[pe_load_hash] %d Original sequences loaded: %.2f sec...\n",
			n_seqs, (float) (clock() - t) / CLOCKS_PER_SEC);

	// Load the hash table itself
	sprintf(hash_fn, "%s.hash", fa_fn);
	opt = (hash_opt*) malloc(sizeof(hash_opt));
	fp = xopen(hash_fn, "rb");
	if (!fread(opt, sizeof(hash_opt), 1, fp)) {
		err_fatal(__func__, "Unable to read from the hash file %s! \n", hash_fn);
	}
	show_msg(
			__func__,
			"Hashing options: k=%d, read_len=%d, n_k_mers=%" ID64 ", n_pos=%" ID64 "...\n",
			opt->k, opt->read_len, opt->n_k_mers, opt->n_pos);
	h->k_mers_occ_acc = (hash_key*) calloc(opt->n_k_mers, sizeof(hash_key));
	h->pos = (hash_value*) calloc(opt->n_pos, sizeof(hash_value));
	h->o = opt;
	fread(h->k_mers_occ_acc, sizeof(hash_key), opt->n_k_mers, fp);
	fread(h->pos, sizeof(hash_value), opt->n_pos, fp);
	fclose(fp);
	show_msg(
			__func__,
			"Hash table loaded, k-mer records: %" ID64 ", positions: %" ID64 " %.2f sec\n",
			opt->n_k_mers, opt->n_pos, (float) (clock() - t) / CLOCKS_PER_SEC);
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

	if (optind + 1 > argc || opt->block_size < 2) {
		return hash_usage();
	}

	pe_hash_core(argv[optind], opt);
	//	pe_hash_test(argv[optind], opt);
	fprintf(stderr, "[pe_hash] Hashing done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
