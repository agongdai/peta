#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "utils.h"

#ifdef HAVE_PTHREAD
#define THREAD_BLOCK_SIZE 1024
#include <pthread.h>
#endif

bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa) {
	bwa_seqio_t *ks;
	ks = bwa_seq_open(fn_fa);
	return ks;
}
