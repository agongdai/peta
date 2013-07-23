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
#include "utils.h"

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

void p_hash_table(hash_table *ht) {
	hash_opt *o = ht->o;
	uint64_t i = 0, j = 0, n_kmers = (1 << (o->k * 2)) + 1, next = 0;
	bwa_seq_t *r = NULL, *kmer = NULL;
	hash_value value = 0;
	hash_key start = 0, end = 0, key = 0;
	int read_id = 0, locus = 0;
	show_debug_msg(__func__, "n_hash_block: %d\n", o->n_hash_block);
	show_debug_msg(__func__, "k: %d\n", o->k);
	show_debug_msg(__func__, "mode: %d\n", o->mode);
	show_debug_msg(__func__, "read_len: %d\n", o->read_len);
	show_debug_msg(__func__, "interleaving: %d\n", o->interleaving);
	show_debug_msg(__func__, "block_size: %d\n", o->block_size);
	for (i = 0; i < n_kmers - 1; i++) {
		next = i + 1;
		start = ht->k_mers_occ_acc[i];
		end = ht->k_mers_occ_acc[next];
		if (end > start) {
			show_debug_msg(__func__, "Kmer: %" ID64 "\n", i);
			kmer = get_kmer_seq(i, o->k);
			p_query("KMER", kmer);
			for (j = start; j < end; j++) {
				value = ht->pos[j];
				read_hash_value(&read_id, &locus, value);
				show_debug_msg(__func__, "Read id: %d; locus: %d\n", read_id,
						locus);
				r = &ht->seqs[read_id];
				r->pos = locus;
				p_query("HIT", r);
			}
		}
	}
}

hash_opt *init_hash_opt() {
	hash_opt *o = (hash_opt*) calloc(1, sizeof(hash_opt));
	o->mode = (BWA_MODE_GAPE | BWA_MODE_COMPREAD);
	o->k = 11;
	o->read_len = 0;
	o->n_k_mers = 0;
	o->n_pos = 0;
	o->interleaving = 2;
	o->n_hash_block = 5;
	o->block_size = o->k / 2;
	return o;
}

uint64_t get_kmer_int(const ubyte_t *seq, const int start,
		const int interleaving, const int len) {
	uint64_t key = 0;
	int i = 0;
	for (i = 0; i < len * interleaving; i += interleaving) {
		if (seq[i + start] < 0 || seq[i + start] > 4) {
			return 0;
		}
		key <<= 2;
		// printf("%d \n", seq[i + start]);
		// Keep the last two digits and perform 'or' oepration
		key = key | (3 & seq[i + start]);
		// printf("%d: key = %" ID64 "\n", i, key);
	}
	return key;
}

hash_key get_hash_key(ubyte_t *seq, const int start, const int interleaving,
		const int len) {
	hash_key key = 0ULL;
	int i = 0;
	for (i = 0; i < len * interleaving; i += interleaving) {
		if (seq[i + start] < 0 || seq[i + start] > 3) {
			seq[i + start] = 0;
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
	//printf("[get_hash_value] (%" ID64 ", %d) => %" ID64 "\n", seq_id, pos_start, value);
	return value;
}

void read_hash_value(index64 *seq_id, int *pos_start, hash_value value) {
	*seq_id = value >> N_POS_BITS;
	*seq_id = *seq_id & HASH_VALUE_HIGHER;
	*pos_start = value & HASH_VALUE_LOWER;
}

void k_hash_core(const char *fa_fn, hash_opt *opt) {
	bwa_seq_t *s, *part_seqs;
	bwa_seqio_t *ks;
	index64 n_k_mers = 0;
	hash_value value;
	int hash_len = 0, tmp = 0, tmp_2 = 0, hash_start = 0, block_no = 0;
	index64 i = 0, i_acc = 0, n_pos = 0, pos_index = 0;
	uint64_t n_part_seqs = 0, n_seqs = 0;
	clock_t t = clock();
	FILE *hash_fp;
	hash_key *k_mers_occ_acc, key = 0ULL;
	hash_value *pos;
	char *hash_fn = malloc(BUFSIZE);
	int *n_occ;

	// All k-mer combinations
	n_k_mers = (1 << (opt->k * 2)) + 1;
	fprintf(stderr, "Hashing library %s ...\n", fa_fn);
	fprintf(stderr,
			"[pe_hash_core] K-mer size: %d; k-mer hashtable size: %" ID64 "\n",
			opt->k, n_k_mers);
	k_mers_occ_acc = (hash_value*) calloc(n_k_mers, sizeof(hash_value));
	n_occ = (int*) calloc(n_k_mers, sizeof(int)); // Store how many occs through the iteration
	hash_len = (opt->read_len / opt->k) * opt->k;

	ks = bwa_open_reads(opt->mode, fa_fn);
	// Round 1: count occurrences of all k-mers
	fprintf(stderr,
			"[pe_hash_core] Round 1/2: Counting occurrences of k-mers ... \n");
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
				k_mers_occ_acc[key]++;
				n_pos++;

				/**
				 printf("key = %" ID64 "\n", key);
				 printf("Hash start: %d\n", hash_start);
				 printf("=== Key %" ID64 ": %" ID64 "===\n", key, k_mers_occ_acc[key]);
				 **/
				key = get_hash_key(s->seq, hash_start + opt->interleaving - 1,
						opt->interleaving, opt->k);
				k_mers_occ_acc[key]++;
				n_pos++;
				/**
				 printf("key = %" ID64 "\n", key);
				 printf("Hash start: %d\n", hash_start);
				 printf("=== Key %" ID64 ": %" ID64 "===\n", key, k_mers_occ_acc[key]);
				 **/
				block_no++;
			}
			hash_start = 0;
			block_no = 0;
		}
		fprintf(stderr,
				"[pe_hash_core] %"ID64" sequences counted: %.2f sec ... \n",
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
			"[pe_hash_core] Pos array size: %" ID64 ", start hashing ...\n",
			n_pos);
	pos = (hash_value*) calloc(n_pos, sizeof(hash_value));

	hash_start = 0;
	block_no = 0;
	fprintf(stderr, "[pe_hash_core] Round 2/2: Store k-mer pointers %s ... \n",
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
	free(hash_fn);
}

hash_table *load_k_hash(char *fa_fn) {
	hash_table *h;
	FILE *fp;
	hash_opt *opt;
	uint32_t n_seqs = 0;
	bwa_seq_t *seqs = NULL;
	char *hash_fn = NULL;
	clock_t t = NULL;
	seqs = load_reads(fa_fn, &n_seqs);
	hash_fn = malloc(FNLEN);
	t = clock();

	fprintf(stderr, "[pe_load_hash] Loading hash table of %s ... \n", fa_fn);
	h = (hash_table*) malloc(sizeof(hash_table));
	h->seqs = seqs;
	h->n_seqs = n_seqs;
	fprintf(stderr,
			"[pe_load_hash] %d Original sequences loaded: %.2f sec ...\n",
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

GPtrArray *find_reads_on_ht(hash_table *ht, bwa_seq_t *query, GPtrArray *hits,
		const int mismatches, const int rev) {
	hash_key key = 0;
	hash_value value = 0;
	int i = 0, j = 0, locus = 0, start = 0, end = 0;
	index64 seq_id = 0;
	bwa_seq_t *r = NULL, *seqs = ht->seqs, *tmp = NULL;
	hash_opt *opt = ht->o;
	int block_no = 0, hash_start = 0;
	if (!hits)
		hits = g_ptr_array_sized_new(0);
	//p_shift_query(query, 0);
	while (block_no < opt->n_hash_block && hash_start <= opt->read_len - opt->k
			* (opt->interleaving)) {
		hash_start = block_no * opt->block_size;
		// Hash two keys for the starting region of a read, interleaving by 1 by default.
		key = get_hash_key(query->seq, hash_start, opt->interleaving, opt->k);

		/**
		 tmp = get_kmer_seq(key, 25);
		 p_query("KMER", tmp);
		 bwa_free_read_seq(1, tmp);
		 **/

		start = ht->k_mers_occ_acc[key];
		end = (key >= opt->n_k_mers) ? ht->k_mers_occ_acc[opt->n_k_mers - 1]
				: ht->k_mers_occ_acc[key + 1];
		if (end > start) {
			for (i = start; i < end; i++) {
				value = ht->pos[i];
				read_hash_value(&seq_id, &locus, value);
				r = &seqs[seq_id];
				if (locus == hash_start) {
					r->pos = (r->pos <= 0) ? 1 : r->pos + 1;
					g_ptr_array_add(hits, r);
				}
				//r->pos = locus;
				//p_shift_query(r, locus);
				//g_ptr_array_add(hits, r);
			}
		}
		block_no++;
	}
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		if (r->pos >= block_no / 2 && head_tail_similar(r, query, opt->k,
				mismatches)) {
			//p_shift_query(r, 0);
			// It is no use for grouping reads;
			// For assembling, it helps because one read is likely to be on one template
			if (rev)
				r->rev_com = 1;
		} else {
			g_ptr_array_remove_index_fast(hits, i--);
		}
		r->pos = -1;
	}
	//show_debug_msg(__func__, "----------------------\n");
	return hits;
}

GPtrArray *find_both_fr_full_reads(hash_table *ht, bwa_seq_t *query,
		GPtrArray *hits, const int mismatches) {
	find_reads_on_ht(ht, query, hits, mismatches, 0);
	set_rev_com(query);
	find_reads_on_ht(ht, query, hits, mismatches, 1);
	return hits;
}

/**
 * Find all reads containing the kmers in seq
 */
GPtrArray *find_reads_with_kmer(hash_table *ht, GPtrArray *hits, int8_t status,
		ubyte_t *seq, index64 len) {
	int64_t i = 0, start = 0, end = 0, seq_id = 0;
	int64_t abs_locus = 0;
	int locus = 0;
	hash_key key = 0;
	hash_value value = 0;
	hash_opt *opt = ht->o;
	bwa_seq_t *r = NULL, *seqs = ht->seqs;
	if (!hits)
		hits = g_ptr_array_sized_new(0);
	// For every possible kmer
	for (i = 0; i <= len - opt->k * opt->interleaving; i++) {
		key = get_hash_key(seq, i, opt->interleaving, opt->k);
		start = ht->k_mers_occ_acc[key];
		end = (key >= opt->n_k_mers) ? ht->k_mers_occ_acc[opt->n_k_mers - 1]
				: ht->k_mers_occ_acc[key + 1];
		if (end > start) {
			for (i = start; i < end; i++) {
				value = ht->pos[i];
				read_hash_value(&seq_id, &locus, value);
				r = &seqs[seq_id];
				abs_locus = locus - i;
				if (abs_locus >= 0 && abs_locus < r->len) {
					if (status == ANY_STATUS || r->status == status)
						g_ptr_array_add(hits, r);
				}
			}
		}
	}
	return hits;
}

/**
 * Align query to the hash_table;
 * The query length is shorter than read length
 */
GPtrArray *align_query(hash_table *ht, GPtrArray *hits, bwa_seq_t *query,
		int8_t status, int mismatches) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	if (!hits)
		hits = g_ptr_array_sized_new(0);
	find_reads_with_kmer(ht, hits, status, query->seq, query->len);
	set_rev_com(query);
	find_reads_with_kmer(ht, hits, status, query->rseq, query->len);
	rm_duplicate(hits);
	return hits;
}

int test_k_hash(char *fa, hash_opt *opt) {
	//k_hash_core(fa, opt);
	hash_table *ht = load_k_hash(fa);
	//p_hash_table(ht);
	bwa_seq_t *query = NULL, *h = NULL;
	GPtrArray *hits = NULL;
	int i = 0, j = 0;
	for (i = 0; i < ht->n_seqs; i++) {
		query = &ht->seqs[i];
		hits = find_reads_on_ht(ht, query, hits, N_MISMATCHES);
		if (hits->len > 100)
			break;
	}
	g_ptr_array_free(hits, TRUE);
}

int k_hash(int argc, char *argv[]) {
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

	k_hash_core(argv[optind], opt);
	//test_k_hash(argv[optind], opt);
	fprintf(stderr, "[pe_hash] Hashing done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
