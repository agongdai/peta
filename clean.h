#include <stdint.h>
#include <inttypes.h>

#define MAX_N16	4294967295
#define N_CHUNK_SEQS	4194304 // # of reads read in every time
#define N_DF_MAX_SEQS	4194304

typedef struct {
	int kmer;
	int mode;
	char *lib_name;
} clean_opt;

typedef struct {
	int key;
	uint16_t count;
} kmer;

typedef struct {
	uint32_t read_id;
	double k_freq;
	double k_sd;
} counter;

int clean_reads(int argc, char *argv[]);
