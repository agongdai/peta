#ifndef K_HASH_H_
#define K_HASH_H_

#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwtaln.h"
#include "utils.h"

#define ID64                  	PRId64
#define N_POS_BITS            	16
#define HASH_VALUE_HIGHER     	281474976710655 	// 48 1's
#define HASH_VALUE_LOWER      	65535				// 16 1's
#define TPL_LOCUS_LOWER			268435455			// 28 1's
#define N_CHUNK_SEQS			4194304 			// # of reads read in every time
#define BWA_MODE				3
#define ANY_STATUS				127					// 8 1's

typedef uint64_t hash_key;
typedef uint64_t hash_value;
typedef int16_t read_pos;

typedef struct {
	int k;				// k value, 25 by default
	int mode; 			// For sequences reading using BWA
	int read_len;		// Read length
	int interleaving; 	// If 2, means hash pattern "10101010...", which is k-weight
	int n_hash_block; 	// How many blocks to hash. Each block varies the size of k/2 from the previous one
	int block_size; 	// How many bases in one block
	index64 n_k_mers;	// # of kmers
	index64 n_pos;		// # of pos list
} hash_opt;

typedef struct {
	hash_opt *o;
	hash_key *k_mers_occ_acc;
	hash_value *pos;
	bwa_seq_t *seqs;
	uint32_t *n_kmers;
	index64 n_seqs;
} hash_table;

#ifdef __cplusplus
extern "C" {
#endif

	int has_next_bit(hash_table *ht, bwa_seq_t *query, int ori);
	hash_table *load_k_hash(char *hash_fn);
	void destroy_ht(hash_table *ht);
	uint64_t get_kmer_int(const ubyte_t *seq, const int start,
			const int interleaving, const int len);
	bwa_seq_t *get_key_seq(uint64_t kmer, const int k);
	hash_value get_hash_value(const index64 seq_id, const int pos_start);
	void read_hash_value(index64 *seq_id, int *pos_start, hash_value value);
	void shrink_ht(hash_table *ht);
	void reload_table(hash_table *ht, char *fa_fn);
	void k_hash_core(const char *fa_fn, hash_opt *opt);
	int k_hash(int argc, char *argv[]);
	hash_key get_hash_key(ubyte_t *seq, const int start,
			const int interleaving, const int len);
	GPtrArray *find_reads_on_ht(hash_table *ht, bwa_seq_t *query, GPtrArray *hits,
			const int mismatches);
	GPtrArray *find_both_fr_full_reads(hash_table *ht, bwa_seq_t *query, GPtrArray *hits,
			const int mismatches);
	GPtrArray *align_query(hash_table *ht, bwa_seq_t *query,
			int8_t status, int mismatches);

#ifdef __cplusplus
}
#endif

#endif
