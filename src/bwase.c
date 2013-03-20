#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stdaln.h"
#include "bwase.h"
#include "bwtaln.h"
#include "bntseq.h"
#include "utils.h"
#include "kstring.h"
#include "pealn.h"

int g_log_n[256];
char *bwa_rg_line, *bwa_rg_id;

void bwa_print_sam_PG();

