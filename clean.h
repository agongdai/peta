#include <stdint.h>

#define ID64 				PRId64

typedef struct {
	int kmer;
	int mode;
	char *lib_name;
} clean_opt;

int clean_reads(int argc, char *argv[]);
