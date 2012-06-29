#include <stdio.h>
#include <string.h>
#include "main.h"
#include "pehash.h"
#include "ass.h"
#include "pealn.h"
#include "eva.h"
#include "clean.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1"
#endif

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: peta (paired-end transcriptome assembler)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Shaojiang Cai <caishaojiang@gmail.com>\n\n");
	fprintf(stderr, "Usage:   peta <command>	[options]\n\n");
	fprintf(stderr, "Command: simu		simulate RNA-seq and RNA-pet sequences\n");
	fprintf(stderr, "Command: aln		align a query to libraries\n");
	fprintf(stderr, "Command: index		create indexes using BWA\n");
	fprintf(stderr, "Command: hash		create a hash table\n");
	fprintf(stderr, "Command: ass		assemble transcriptomes\n");
	fprintf(stderr, "Command: eva		evaluate the performance\n");
	fprintf(stderr, "\n");
	return 1;
}

void bwa_print_sam_PG() {
	//	printf("@PG\tID:peta\tPN:peta\tVN:%s\n", PACKAGE_VERSION);
}

int test() {

	return 0;
}

int main(int argc, char *argv[]) {
	if (argc < 2)
		return usage();
	if (strcmp(argv[1], "simu") == 0)
		return peta_simu(argc - 1, argv + 1);
	else if (strcmp(argv[1], "index") == 0)
		return bwa_index(argc - 1, argv + 1);
	else if (strcmp(argv[1], "hash") == 0)
		return pe_hash(argc - 1, argv + 1);
	else if (strcmp(argv[1], "bwa_aln") == 0)
		return bwa_aln(argc - 1, argv + 1);
	else if (strcmp(argv[1], "pe_aln") == 0)
			return pe_aln_test(argc - 1, argv + 1);
	else if (strcmp(argv[1], "ass") == 0)
		return pe_ass(argc - 1, argv + 1);
	else if (strcmp(argv[1], "samse") == 0)
		return bwa_sai2sam_se(argc - 1, argv + 1);
	else if (strcmp(argv[1], "eva") == 0)
		return eva_main(argc - 1, argv + 1);
	else if (strcmp(argv[1], "test") == 0)
		return test(argc - 1, argv + 1);
	else if (strcmp(argv[1], "clean") == 0)
		return clean_reads(argc - 1, argv + 1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
