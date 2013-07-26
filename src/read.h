#ifndef READ_H_
#define READ_H_

#include <stdint.h>
#include <inttypes.h>
#include <glib.h>
#include "bwtaln.h"
#include "utils.h"
#include "k_hash.h"

typedef struct {
	uint32_t *similar_reads_count;
	index64 n_seqs;
} read_hash;

typedef struct {
	hash_table *ht;
	read_hash *rh;
} align_para;

#ifdef __cplusplus
extern "C" {
#endif

	void destroy_rh(read_hash *rh);
	int group_main(int argc, char *argv[]);
	read_hash *load_read_hash(char *fa);

#ifdef __cplusplus
}
#endif

#endif
