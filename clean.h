#include <stdint.h>

typedef struct {
	int kmer;
	int mode;
	char *lib_name;
} clean_opt;

typedef struct {
	int key;
	int count;
} kmer;

int clean_reads(int argc, char *argv[]);
