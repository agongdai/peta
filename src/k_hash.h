#ifndef K_HASH_H_
#define K_HASH_H_

#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwase.h"
#include "utils.h"

#define ID64                    PRId64
#define K_N_POS_BITS            64
#define K_HASH_VALUE_HIGHER     281474976710655 // 48 1's
#define K_HASH_VALUE_LOWER      65535

typedef uint64_t hash_key;
typedef uint64_t hash_value;
typedef int16_t read_pos;

typedef struct {
    int k;
    int mode;           // For sequences reading using BWA
    int read_len;
    int interleaving;   // If 2, means hash pattern "10101010...", which is k-weight
    int n_hash_block;   // How many blocks to hash. Each block varies the size of k/2 from the previous one
    int block_size;     // How many bases in one block
    uint64_t n_k_mers;
    uint64_t n_pos;
} hash_opt;

typedef struct {
    hash_opt *o;
    hash_key *k_mers_occ_acc;
    hash_value *pos;
    bwa_seq_t *seqs;
    uint64_t n_seqs;
} hash_table;

hash_table *load_k_hash(char *hash_fn);
