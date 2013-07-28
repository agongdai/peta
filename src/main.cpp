/*

 * kmer.cpp
 *
 *  Created on: 02-Apr-2013
 *      Author: carl
 */

#include <stdio.h>
#include <string.h>
#include <glib.h>
#include "main.h"
#include "ass.hpp"
#include "kmers.hpp"
#include "oracle.hpp"
#include "k_hash.h"

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
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[]) {
	//	test_kmer_hash(
	//			"/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897_corrected.fa");
	//	build_kmers_hash("/home/carl/Projects/peta/rnaseq/hg19/SRX011545/SRR027876.fa", 25, 1);
	//	return 1;
	if (!g_thread_supported())
		g_thread_init( NULL);
	if (argc < 2)
		return usage();
	else if (strcmp(argv[1], "ass") == 0)
		return pe_kmer(argc - 1, argv + 1);
	else if (strcmp(argv[1], "hash") == 0)
		return build_kmer_hash(argc - 1, argv + 1);
	else if (strcmp(argv[1], "oracle") == 0)
		return oracle_set(argc - 1, argv + 1);
	else if (strcmp(argv[1], "exon") == 0)
		return genome_splicings(argc - 1, argv + 1);
	else if (strcmp(argv[1], "freq") == 0)
		return export_frequency(argc - 1, argv + 1);
	else if (strcmp(argv[1], "k_hash") == 0)
		return k_hash(argc - 1, argv + 1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}

