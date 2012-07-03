#include <stdint.h>
#include <inttypes.h>

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
	char *read_id;
	int k_freq;
} counter;

int clean_reads(int argc, char *argv[]);
