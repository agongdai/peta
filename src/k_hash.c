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

int n_threads = 4;

/**
 * Check whether the query has next character
 */
int has_next_bit(hash_table *ht, bwa_seq_t *query, int ori) {
	uint64_t key = 0, new_key = 0;
	int i = 0, start = 0, end = 0, n_multi = 0;
	hash_opt *opt = ht->o;
	key = get_hash_key(query->seq, i, opt->interleaving, opt->k);
	for (i = 0; i < 4; i++) {
		new_key = shift_bit(key, i, ht->o->k, ori);

		start = ht->k_mers_occ_acc[new_key];
		end
				= (new_key >= opt->n_k_mers) ? ht->k_mers_occ_acc[opt->n_k_mers
						- 1] : ht->k_mers_occ_acc[new_key + 1];
		if (end > start + 1) {
			n_multi++;
		}
		/**
		 show_debug_msg(__func__, "---\n");
		 bwa_seq_t *key_seq = get_key_seq(key, 11);
		 p_query(__func__, key_seq);
		 bwa_free_read_seq(1, key_seq);
		 show_debug_msg(__func__, "Next: %c; shift to %s; Multi: %d\n",
		 "ACGTN"[i], ori ? "left" : "right", n_multi);
		 key_seq = get_key_seq(new_key, 11);
		 p_query(__func__, key_seq);
		 bwa_free_read_seq(1, key_seq);
		 **/
		if (n_multi >= 2)
			return 1;
	}
	return 0;
}

bwa_seq_t *get_key_seq(uint64_t kmer, const int k) {
	ubyte_t *seq = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	uint64_t copy = kmer, index = 0;
	read = blank_seq(k);
	seq = read->seq;
	for (i = 0; i < k; i++) {
		index = copy;
		index &= 3; // Keep only the last two bits
		seq[k - 1 - i] = index;
		copy >>= 2;
	}
	read->len = k;
	read->status = FRESH;
	read->name = (char*) malloc(64);
	sprintf(read->name, "%" ID64, kmer);
	set_rev_com(read);
	return read;
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
		free(ht->n_kmers);
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
			kmer = get_key_seq(i, o->k);
			p_query("KMER", kmer);
			for (j = start; j < end; j++) {
				value = ht->pos[j];
				read_hash_value(&read_id, &locus, value);
				show_debug_msg(__func__, "Hit id: %d; locus: %d\n", read_id, locus);
				if (ht->seqs) {
					r = &ht->seqs[read_id];
					r->pos = locus;
					p_query("HIT", r);
				}
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


/**
 * Hash all k-mers on a template
 * Assumption: the memory hold by *hash must be at least 4^k * sizeof(uint32_t)
 */
void hash_a_tpl(uint32_t *hash, const ubyte_t *seq, const int len, const int k, const int tpl_id) {
	uint32_t i = 0;
	hash_key key = 0ULL;
	hash_value value = 0ULL;
	if (len < k || !hash || tpl_id < 0 || !seq || k <= 0) return;
	for (i = 0; i <= len - k; i++) {
		key = get_hash_key(seq[i], i, 1, k);
		value = get_hash_value(tpl_id, i);
		hash[key] = value;
	}
}

void read_hash_value(index64 *seq_id, int *pos_start, hash_value value) {
	*seq_id = value >> N_POS_BITS;
	*seq_id = *seq_id & HASH_VALUE_HIGHER;
	*pos_start = value & HASH_VALUE_LOWER;
}

/**
 * Remove the USED/DEAD read occurrences from hash table
 */
void shrink_ht(hash_table *ht) {
	int i = 0, j = 0, n_valid = 0, n_empty = 0, pos_index = 0;
	index64 n_k_mers = 0, seq_id = 0, locus = 0;
	hash_opt *opt = NULL;
	hash_key *k_mers_occ_acc = NULL, *new_k_mers = NULL, key_this = 0,
			key_next = 0;
	hash_value *pos = NULL, *new_pos = NULL, value = 0;
	bwa_seq_t *r = NULL;

	opt = ht->o;
	k_mers_occ_acc = ht->k_mers_occ_acc;
	pos = ht->pos;
	n_k_mers = (1 << (opt->k * 2)) + 1;
	new_k_mers = (hash_value*) calloc(n_k_mers, sizeof(hash_value));
	for (i = 0; i < n_k_mers - 1; i++) {
		key_this = k_mers_occ_acc[i];
		key_next = k_mers_occ_acc[i + 1];
		for (j = key_this; j < key_next; j++) {
			value = pos[j];
			read_hash_value(&seq_id, &locus, value);
			if (seq_id >= ht->n_seqs || seq_id < 0) {
				show_debug_msg(__func__,
						"WARNING: Sequence id not correct: %d \n", seq_id);
				continue;
			}
			r = &ht->seqs[seq_id];
			if (value == 0 || r->status == USED || r->status == DEAD) {
				n_empty++;
			} else {
				n_valid++;
			}
		}
		new_k_mers[i + 1] = key_next - n_empty;
	}
	free(k_mers_occ_acc);
	ht->k_mers_occ_acc = new_k_mers;

	new_pos = (hash_value*) calloc(n_valid + 1, sizeof(hash_value));
	for (i = 0; i < n_valid + n_empty; i++) {
		value = pos[i];
		read_hash_value(&seq_id, &locus, value);
		if (seq_id >= ht->n_seqs || seq_id < 0) {
			show_debug_msg(__func__, "WARNING: Sequence id not correct: %d \n",
					seq_id);
			continue;
		}
		r = &ht->seqs[seq_id];
		if (value > 0 && r->status != USED && r->status != DEAD) {
			new_pos[pos_index] = value;
			pos_index++;
		}
	}
	free(pos);
	ht->pos = new_pos;
	opt->n_pos = n_valid;
}

void k_hash_core(const char *fa_fn, hash_opt *opt) {
	bwa_seq_t *s, *part_seqs;
	bwa_seq_t *key_seq = NULL;
	bwa_seqio_t *ks;
	index64 n_k_mers = 0;
	hash_value value;
	int hash_len = 0, tmp = 0, tmp_2 = 0, hash_start = 0, block_no = 0;
	index64 i = 0, i_acc = 0, n_pos = 0, pos_index = 0;
	uint64_t n_part_seqs = 0, n_seqs = 0;
	clock_t t = clock();
	FILE *hash_fp = NULL;
	hash_key *k_mers_occ_acc, key = 0ULL;
	hash_value *pos;
	char *hash_fn = (char*) malloc(BUFSIZE);
	int *n_occ;
	uint32_t *kmer_occ_on_reads = NULL;
	GThreadPool *thread_pool = NULL;

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

			// Ignore the reads all 'AAAAATAAAA', etc
			//if (same_bytes(s->seq, s->len) || too_many_ns(s->seq, s->len))
			//	continue;

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
	kmer_occ_on_reads = (uint32_t*) calloc(n_seqs, sizeof(uint32_t));

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

			//if (same_bytes(s->seq, s->len) || too_many_ns(s->seq, s->len))
			//	continue;

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
				"[pe_hash_core] %"ID64" sequences hashed: %.2f sec ... \n",
				n_seqs, (float) (clock() - t) / CLOCKS_PER_SEC);
		bwa_free_read_seq(n_part_seqs, part_seqs);
	}
	bwa_seq_close(ks);

	sprintf(hash_fn, "%s.hash", fa_fn);
	hash_fp = xopen(hash_fn, "w");
	opt->n_k_mers = n_k_mers;
	opt->n_pos = n_pos;
	fprintf(stderr, "[pe_hash_core] Saving to %s ... \n", hash_fn);
	fwrite(opt, sizeof(hash_opt), 1, hash_fp);
	fprintf(stderr, "[pe_hash_core] Saving k-mers ... \n");
	fwrite(k_mers_occ_acc, sizeof(hash_key), n_k_mers, hash_fp);
	fprintf(stderr, "[pe_hash_core] Saving occurrences ... \n");
	fwrite(pos, sizeof(hash_value), n_pos, hash_fp);

	fprintf(stderr, "[pe_hash_core] Saving kmer counts on reads ... \n");
	// In first round, just store all zero values
	fwrite(kmer_occ_on_reads, sizeof(uint32_t), n_seqs, hash_fp);
	free(kmer_occ_on_reads);

	fclose(hash_fp);
	free(k_mers_occ_acc);
	free(n_occ);
	free(hash_fn);
}

hash_table *load_k_hash(char *fa_fn) {
	hash_table *ht;
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
	ht = (hash_table*) malloc(sizeof(hash_table));
	ht->seqs = seqs;
	ht->n_seqs = n_seqs;
	ht->n_kmers = (uint32_t*) calloc(n_seqs, sizeof(uint32_t));
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
	show_msg(__func__,
			"Hashing options: k=%d, read_len=%d, n_k_mers=%" ID64 ", n_pos=%" ID64 " ...\n",
			opt->k, opt->read_len, opt->n_k_mers, opt->n_pos);
	ht->k_mers_occ_acc = (hash_key*) calloc(opt->n_k_mers, sizeof(hash_key));
	ht->pos = (hash_value*) calloc(opt->n_pos, sizeof(hash_value));
	ht->o = opt;
	fread(ht->k_mers_occ_acc, sizeof(hash_key), opt->n_k_mers, fp);
	fread(ht->pos, sizeof(hash_value), opt->n_pos, fp);
	fread(ht->n_kmers, sizeof(uint32_t), n_seqs, fp);
	fclose(fp);
	show_msg(__func__,
			"Hash table loaded, k-mer records: %" ID64 ", positions: %" ID64 " %.2f sec\n",
			opt->n_k_mers, opt->n_pos, (float) (clock() - t) / CLOCKS_PER_SEC);
	free(hash_fn);
	return ht;
}

/**
 * Reload the hashtable except the sequences
 */
void reload_table(hash_table *ht, char *fa_fn) {
	char *hash_fn = malloc(FNLEN);
	sprintf(hash_fn, "%s.hash", fa_fn);
	FILE *fp = xopen(hash_fn, "rb");
	free(ht->o);
	ht->o = (hash_opt*) malloc(sizeof(hash_opt));
	if (!fread(ht->o, sizeof(hash_opt), 1, fp)) {
		err_fatal(__func__, "Unable to read from the hash file %s! \n", hash_fn);
	}
	free(ht->k_mers_occ_acc);
	free(ht->pos);
	ht->k_mers_occ_acc = (hash_key*) calloc(ht->o->n_k_mers, sizeof(hash_key));
	ht->pos = (hash_value*) calloc(ht->o->n_pos, sizeof(hash_value));
	fread(ht->k_mers_occ_acc, sizeof(hash_key), ht->o->n_k_mers, fp);
	fread(ht->pos, sizeof(hash_value), ht->o->n_pos, fp);
	free(hash_fn);
	fclose(fp);
}

/**
 * The query length is the same as read length.
 * Return a read half of the hashed kmers are matched.
 */
GPtrArray *find_reads_on_ht(hash_table *ht, bwa_seq_t *query, GPtrArray *hits,
		const int mismatches) {
	hash_key key = 0;
	hash_value value = 0;
	int i = 0, j = 0, locus = 0, start = 0, end = 0, n_block_kmers = 2;
	index64 seq_id = 0;
	bwa_seq_t *r = NULL, *seqs = ht->seqs, *tmp = NULL;
	hash_opt *opt = ht->o;
	int block_no = 0, hash_start = 0;
	if (!hits)
		hits = g_ptr_array_sized_new(0);
	//p_shift_query(query, 0);
	while (block_no < opt->n_hash_block && hash_start <= opt->read_len - opt->k
			* (opt->interleaving)) {
		// Hash two keys for the starting region of a read, interleaving by 1 by default.
		for (j = 0; j < n_block_kmers; j++) {
			hash_start = block_no * opt->block_size + j;
			if (hash_start <= opt->read_len - opt->k * (opt->interleaving)) {
				key = get_hash_key(query->seq, hash_start, opt->interleaving,
						opt->k);

				start = ht->k_mers_occ_acc[key];
				end = (key >= opt->n_k_mers) ? ht->k_mers_occ_acc[opt->n_k_mers
						- 1] : ht->k_mers_occ_acc[key + 1];
				for (i = start; i < end; i++) {
					value = ht->pos[i];
					read_hash_value(&seq_id, &locus, value);
					r = &seqs[seq_id];
					//if (strcmp(r->name, "15398") == 0) {
					//	p_query(__func__, r);
					//}
					if (locus == hash_start) {
						// To avoid too many duplicates in hits (some extreme case)
						if (r->pos == IMPOSSIBLE_NEGATIVE)
							g_ptr_array_add(hits, r);
						//p_query(__func__, r);
						// Here 'pos' stores how many kmers this query contains
						r->pos = (r->pos <= 0) ? 1 : r->pos + 1;
					}
				}
			}
		}
		block_no++;
	}
	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		// p_query(__func__, r);
		// If half of the kmers are obtained, go for head_tail_similar in find_both_fr_full_reads
		if (r->pos < block_no) {
			// Remove the read from hits and reset the pos
			g_ptr_array_remove_index_fast(hits, i--);
			r->pos = IMPOSSIBLE_NEGATIVE;
		}
	}
	//show_debug_msg(__func__, "----------------------\n");
	return hits;
}

/**
 * The query length is exactly the same as read length
 * Report a hit if a read:
 * 	1. Half of the hashed kmers are present;
 *  2. Head and tail 11 bases are similar
 */
GPtrArray *find_both_fr_full_reads(hash_table *ht, bwa_seq_t *query,
		GPtrArray *hits, const int mismatches) {
	index64 i = 0;
	int rev_com = 0;
	bwa_seq_t *r = NULL;
	//if (strcmp(query->name, "SOME") == 0)
	//	p_query(__func__, query);
	hits = find_reads_on_ht(ht, query, hits, mismatches);
	switch_fr(query);
	hits = find_reads_on_ht(ht, query, hits, mismatches);
	switch_fr(query);

	for (i = 0; i < hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(hits, i);
		//if (query->full_len == 1000)
			//p_query(__func__, r);
		r->pos = IMPOSSIBLE_NEGATIVE;
		if (head_tail_similar(r, query, ht->o->k, mismatches, &rev_com)) {
			r->rev_com = rev_com;
		} else {
			g_ptr_array_remove_index_fast(hits, i--);
		}
	}
	return hits;
}

/**
 * Find all reads containing the kmers in seq
 */
GPtrArray *find_reads_with_kmer(hash_table *ht, GPtrArray *hits, int8_t status,
		ubyte_t *seq, index64 len) {
	int64_t i = 0, j = 0, start = 0, end = 0, seq_id = 0;
	int64_t abs_locus = 0;
	int locus = 0;
	uint64_t all_ones = 1;
	hash_key even_key = 0, odd_key = 0, key = 0;
	hash_value value = 0;
	hash_opt *opt = ht->o;
	bwa_seq_t *r = NULL, *seqs = ht->seqs;
	// opt->k 1's
	all_ones <<= opt->k * 2;
	all_ones--;
	even_key = get_hash_key(seq, 0, opt->interleaving, opt->k);
	odd_key = get_hash_key(seq, 1, opt->interleaving, opt->k);
	for (i = 0; i <= len - opt->k * opt->interleaving + 1; i++) {
		// Shift the kmer to left
		if (i % 2 == 0) {
			if (i >= 2) {
				even_key <<= 2;
				even_key |= (3 & seq[i + (opt->k - 1) * opt->interleaving]);
				even_key &= all_ones;
			}
			key = even_key;
		} else {
			if (i >= 2) {
				odd_key <<= 2;
				odd_key |= (3 & seq[i + (opt->k - 1) * opt->interleaving]);
				odd_key &= all_ones;
			}
			key = odd_key;
		}
		//show_debug_msg(__func__, "KEY: %" ID64 ". \n", key);
		/**
		 if (seq[0] == 3 && seq[1] == 1 && seq[2] == 3 && seq[3] == 3 && len == 31 && status == FRESH) {
		 show_debug_msg(__func__, "---\n");
		 bwa_seq_t *key_seq = get_key_seq(key, 11);
		 p_query(__func__, key_seq);
		 bwa_free_read_seq(1, key_seq);
		 }**/

		start = ht->k_mers_occ_acc[key];
		end = (key >= opt->n_k_mers) ? ht->k_mers_occ_acc[opt->n_k_mers - 1]
				: ht->k_mers_occ_acc[key + 1];
		//show_debug_msg(__func__, "Start~End: [%" ID64 ", %" ID64 "]\n", start, end);
		for (j = start; j < end; j++) {
			value = ht->pos[j];
			read_hash_value(&seq_id, &locus, value);
			r = &seqs[seq_id];
			// To avoid duplicate adding
			if (r->pos != IMPOSSIBLE_NEGATIVE)
				continue;
			abs_locus = locus - i;
			/**
			 if (seq[0] == 3 && seq[1] == 1 && seq[2] == 3 && seq[3] == 3 && len == 31 && status == FRESH) {
			 p_query(__func__, r);
			 show_debug_msg(__func__, "i: %d; locus: %d; abs_locus: %d \n", i, locus, abs_locus);
			 }**/
			if (abs_locus >= 0 && abs_locus <= r->len + 1 - opt->interleaving
					* opt->k) {
				// If the status of read is as requested
				if (r->status == status || status == ANY_STATUS) {
					//p_query(__func__, r);
					//show_debug_msg(__func__, "Locus: %d\n", locus);
					r->pos = abs_locus;
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
GPtrArray *align_query(hash_table *ht, bwa_seq_t *query, int8_t status,
		int mismatches) {
	index64 i = 0;
	bwa_seq_t *r = NULL;
	GPtrArray *hits = g_ptr_array_sized_new(0), *rev_hits = g_ptr_array_sized_new(0);

	//	if (query->full_len == 1000) {
	//		show_debug_msg(__func__, "Before align \n");
	//		p_query(__func__, query);
	//		bwa_seq_t *rev = new_seq(query, query->len, 0);
	//		switch_fr(rev);
	//		p_query("REV", rev);
	//		bwa_free_read_seq(1, rev);
	//		p_readarray(hits, 1);
	//	}

	hits = find_reads_with_kmer(ht, hits, status, query->seq, query->len);

	//	if (query->full_len == 1000) {
	//		p_query(__func__, query);
	//		bwa_seq_t *rev = new_seq(query, query->len, 0);
	//		switch_fr(rev);
	//		p_query("REV", rev);
	//		bwa_free_read_seq(1, rev);
	//		for (i = 0; i < hits->len; i++) {
	//			r = (bwa_seq_t*) g_ptr_array_index(hits, i);
	//			p_query("SOME", r);
	//		}
	//	}

	set_rev_com(query);
	rev_hits = find_reads_with_kmer(ht, rev_hits, status, query->rseq, query->len);

	// Reset the rev_com first, later it is used as an indicator
	for (i = 0; i < rev_hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(rev_hits, i);
		r->rev_com = 0;
	}

	// If reverse complement, update the pos
	for (i = 0; i < rev_hits->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(rev_hits, i);
		//if (query->full_len == 1000)
		//	p_query(__func__, r);
		// To make sure the 'pos' will be changed once only
		if (!r->rev_com) {
			//if (strcmp(r->name, "15398") == 0)
			//	p_query(__func__, r);
			r->rev_com = 1;
			// Adjust the locus due to reverse complement
			//show_debug_msg(__func__, "Query->len: %d; r->pos: %d \n", query->len, r->pos);
			//p_query("REV_COM", r);
			r->pos = (ht->o->read_len - query->len) - r->pos;
			//p_query("REV_COM", r);
			g_ptr_array_add(hits, r);
		}
	}
	g_ptr_array_free(rev_hits, TRUE);
	hits = rm_duplicates(hits);
	return hits;
}

/**
 * The thread to align one read to the RNA-seq library
 */
gint align_read_thread(gpointer r, gpointer para) {
	bwa_seq_t *query = (bwa_seq_t*) r;
	hash_table *ht = para;
	GPtrArray *hits = g_ptr_array_sized_new(0);
	index64 i = 0;
	bwa_seq_t *tmp = NULL;
	uint32_t *n_kmers = ht->n_kmers;
	index64 read_id = 0;
	read_id = atoll(query->name);
	if (read_id < 0 || read_id >= ht->n_seqs || query->pos
			!= IMPOSSIBLE_NEGATIVE)
		return 0;
	if (is_biased_q(query) || too_many_ns(query, query->len))
		return 0;
	find_both_fr_full_reads(ht, query, hits, N_MISMATCHES);
	for (i = 0; i < hits->len; i++) {
		tmp = (bwa_seq_t*) g_ptr_array_index(hits, i);
		if (seq_ol(tmp, query, tmp->len, N_MISMATCHES) == -1) {
			g_ptr_array_remove_index_fast(hits, i--);
		}
	}
	// Mark 'pos' as not -1, to save time
	//show_debug_msg(__func__, "Hits: %d\n", hits->len);
	for (i = 0; i < hits->len; i++) {
		tmp = (bwa_seq_t*) g_ptr_array_index(hits, i);
		tmp->pos = hits->len;
		n_kmers[atoll(tmp->name)] = hits->len;
	}
	//show_debug_msg(__func__, "---------------------\n");
	n_kmers[read_id] = hits->len;
	g_ptr_array_free(hits, TRUE);
	return 0;
}

void group_reads(hash_table *ht, const int n_threads) {
	GThreadPool *thread_pool = NULL;
	int i = 0;
	bwa_seq_t *r = NULL;

	show_msg(__func__, "Aligning %d reads on %d threads ...\n", ht->n_seqs,
			n_threads);
	thread_pool = g_thread_pool_new((GFunc) align_read_thread, ht, n_threads,
			TRUE, NULL);
	for (i = 0; i < ht->n_seqs; i++) {
		r = (bwa_seq_t*) &ht->seqs[i];
		ht->n_kmers[atoll(r->name)] = 0;
	}
	for (i = 0; i < ht->n_seqs; i++) {
		r = (bwa_seq_t*) &ht->seqs[i];
		//align_read_thread((gpointer)r, (gpointer)para);
		g_thread_pool_push(thread_pool, (gpointer) r, NULL);
		//if (i > 100)
		//	break;
	}
	g_thread_pool_free(thread_pool, 0, 1);
}

void re_hash(const char *fa_fn) {
	char rehash_fn[BUFSIZ];
	FILE *fp = NULL;

	hash_table *ht = load_k_hash(fa_fn);
	group_reads(ht, n_threads);

	sprintf(rehash_fn, "%s.hash", fa_fn);
	fp = xopen(rehash_fn, "w");

	fprintf(stderr, "[re_hash] Saving to %s ... \n", rehash_fn);
	fwrite(ht->o, sizeof(hash_opt), 1, fp);
	fprintf(stderr, "[re_hash] Saving k-mers ... \n");
	fwrite(ht->k_mers_occ_acc, sizeof(hash_key), ht->o->n_k_mers, fp);
	fprintf(stderr, "[re_hash] Saving occurrences ... \n");
	fwrite(ht->pos, sizeof(hash_value), ht->o->n_pos, fp);

	fprintf(stderr, "[re_hash] Saving kmer counts on reads ... \n");
	fwrite(ht->n_kmers, sizeof(uint32_t), ht->n_seqs, fp);
	fclose(fp);
}

int test_k_hash(char *fa, hash_opt *opt) {
	//k_hash_core(fa, opt);
	hash_table *ht = load_k_hash(fa);
	//p_hash_table(ht);
	bwa_seq_t *query = NULL, *h = NULL;
	int i = 0, j = 0;
	for (i = 0; i < ht->n_seqs; i++) {
		query = &ht->seqs[99];
		p_query(__func__, query);
		show_debug_msg(__func__, "HITS: %d \n", ht->n_kmers[99]);
		break;
	}
}

int k_hash(int argc, char *argv[]) {
	int c;
	clock_t t;
	t = clock();
	hash_opt *opt = init_hash_opt();

	while ((c = getopt(argc, argv, "k:l:i:b:s:t")) >= 0) {
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
		case 't':
			n_threads = atoi(optarg);
			break;
		default:
			return 1;
		}
	}

	if (optind + 1 > argc || opt->block_size < 2) {
		return hash_usage();
	}

	if (!g_thread_supported())
		g_thread_init(NULL);

	k_hash_core(argv[optind], opt);
	re_hash(argv[optind]);
	//test_k_hash(argv[optind], opt);
	fprintf(stderr, "[pe_hash] Hashing done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	return 0;
}
