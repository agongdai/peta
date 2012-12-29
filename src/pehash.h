/*
 * pehash.h
 *
 *  Created on: 11-Dec-2011
 *      Author: carl
 */
#ifndef PEHASH_H_
#define PEHASH_H_

#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwase.h"
#include "utils.h"

#define ID64 				PRId64
#define N_POS_BITS			16
#define HASH_VALUE_HIGHER	281474976710655 // 48 1's
#define HASH_VALUE_LOWER	65535

typedef uint64_t hash_key;
typedef uint64_t hash_value;
typedef int16_t read_pos;

typedef struct {
	int k;
	int mode; // For sequences reading using BWA
	int read_len;
	int interleaving; // If 2, means hash pattern "10101010...", which is k-weight
	int n_hash_block; // How many blocks to hash. Each block varies the size of k/2 from the previous one
	int block_size; // How many bases in one block
	index64 n_k_mers;
	index64 n_pos;
} hash_opt;

typedef struct {
	hash_opt *o;
	hash_key *k_mers_occ_acc;
	hash_value *pos;
	bwa_seq_t *seqs;
	index64 n_seqs;
} hash_table;

typedef struct {
	int k;
	index64 n_reads;
	GPtrArray *pos;
} reads_ht;

int pe_hash(int argc, char *argv[]);
hash_table *pe_load_hash(char *hash_fn);
void pe_hash_core(const char *fa_fn, hash_opt *opt);
hash_key get_hash_key(const ubyte_t *seq, const int start,
		const int interleaving, const int len);
hash_value get_hash_value(const index64 seq_id, const int pos_start);
void read_hash_value(index64 *seq_id, int *pos_start, hash_value value);
void destroy_ht(hash_table *ht);
hash_opt *init_hash_opt();

GPtrArray *find_reads_ol_template(reads_ht *ht, bwa_seq_t *template,
		bwa_seq_t *seqs, const int ori);
GPtrArray *find_edges_ol(reads_ht *ht, bwa_seq_t *template,
		GPtrArray *all_edges);
void destroy_reads_ht(reads_ht *ht);
reads_ht *build_reads_ht(const int k, GPtrArray *initial_reads);
reads_ht *build_edges_ht(const int k, GPtrArray *init_edges);
void add_read_to_ht(reads_ht *ht, bwa_seq_t *read);
void rm_read_from_ht(reads_ht *ht, bwa_seq_t *read);

#endif /* PEHASH_H_ */
