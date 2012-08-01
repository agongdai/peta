#include <stdint.h>
#include <inttypes.h>

#define MAX_N16		4294967295
#define LOW_KMER	2
#define UNEVEN_THRE	3

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
	double k_mean;
} counter;

int clean_reads(int argc, char *argv[]);
