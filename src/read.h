#ifndef READ_H_
#define READ_H_

#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwtaln.h"
#include "utils.h"
#include "k_hash.h"

typedef struct {
	GPtrArray *read_groups;
	bwa_seq_t *seqs;
	index64 n_seqs;
} read_hash;

typedef struct {
	hash_table *ht;
	read_hash *rh;
} align_para;

#ifdef __cplusplus
extern "C" {
#endif

	int group_main(int argc, char *argv[]);
	read_hash *load_read_hash(char *fa);

#ifdef __cplusplus
}
#endif

#endif
