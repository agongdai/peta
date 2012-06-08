#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <zlib.h>
#include <glib.h>
#include "main.h"
#include "utils.h"
#include "rand.h"
#include "peseq.h"
#include "pechar.h"

#define GENE_N 60000

float err = 0.02;
int c, mu = 500, sigma = 10, coverage = 50, seq_id = 0, pet_id = 0, tx_id,
		both_files = 0, min_l = 1000000, max_l = 0, ave_l = 0, total_l = 0,
		n_seqs = 0, n_pets = 0, n_tx = 0, n_exons = 0, opt_n50 = 0,
		total_tx_len = 0, read_len = 0, pet_len = 0, strand_specific = 0;
char *rnapet_fn, *tx_fn;
GArray *all_len;

int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   peta simu [-c <coverage>] [-g <gene list>] ");
	fprintf(stderr, "[-err <error rate>] [-m <mu of normal distribution>] ");
	fprintf(stderr, "[-s <sigma of normal distribution>] "
		"[-t <transcript file>] [-p <rnapet file>] <seq.fasta>\n\n");
	return 1;
}

gint cmp_lens(gpointer a, gpointer b) {
	int* x = (int*) a;
	int* y = (int*) b;
	return *y - *x;
}

void rev_cmp(char *seq, int len) {
	int i = 0;
	for (i = 0; i < len >> 1; ++i) {
		char tmp = compl(seq[len - 1 - i]);
		seq[len - 1 - i] = compl(seq[i]);
		seq[i] = tmp;
	}
}

/**
 * Simulate the transcript of one gene
 *
 * tx_fp: file to store the transcript;
 * tx: transcript sequence;
 * fa: original sequence;
 * gn: gene information
 */
void stx(FILE *tx_fp, seq *tx, const seq *fa, gene *gn, int renew) {
	// Get all exon start positions
	int i = 0, j, index_s, index_e, gene_l = 0;
	char *exlist_s[BUFSIZE], *exlist_e[BUFSIZE];
	char *header = malloc(BUFSIZE);
	exlist_s[0] = strtok(gn->exon_s, ",");
	while (exlist_s[i] != NULL) { //ensure a pointer was found
		exlist_s[++i] = strtok(NULL, ","); //continue to tokenize the string
	}

	// Get all exon end positions
	i = 0;
	exlist_e[0] = strtok(gn->exon_e, ",");
	while (exlist_e[i] != NULL) { //ensure a pointer was found
		exlist_e[++i] = strtok(NULL, ","); //continue to tokenize the string
	}

	strcpy(tx->h, gn->id);
	// Generate the transcripts.
	tx->s = (char*) malloc(1024);
	tx->l = 0;
	tx->m = 1024;
	//		tx->h = header;
	//fprintf(stderr, "# of exons: %d\n", gn->exon_c);
	n_exons += gn->exon_c;
	for (i = 0; i < gn->exon_c; i++) {
		index_s = atoi(exlist_s[i]);
		index_e = atoi(exlist_e[i]);
		total_tx_len += index_e - index_s;
		gene_l += index_e - index_s;

		// printf("Exon range of %s: [%d, %d] (%d) \n", attrs[1], index_s,
		// index_e, exon_l);
		for (j = index_s; j < index_e; j++) {
			// Allocate more memory if not sufficient.
			if (tx->l + 1 >= tx->m) {
				tx->m = tx->l + 2;
				kroundup32(tx->m);
				// printf("tx->l = %zd, exon_l = %d, tx->m = %zd\n", tx->l,
				// exon_l, tx->m);
				tx->s = (char*) realloc(tx->s, tx->m);
			}
			if (gn->strand == '+') {
				tx->s[tx->l++] = to_upper_lower(fa->s[j]);
			} else {
				// If "-", get the reverse complement.
				// printf("[index_s, index_e, j] = [%d, %d, %d]\n", index_s, index_e, j);
				// printf("fa->s[index_e + index_s - j] = %c\n", fa->s[index_e + index_s - j]);
				tx->s[tx->l++] = to_upper_lower(compl(fa->s[index_e + index_s
						- j]));
			}
		}
	}
	g_array_append_val(all_len, gene_l);

	// Output the sequence
	if (renew) {
		// Save the header
		sprintf(header, ">%d|%s|len=%zd", tx_id++, gn->id, tx->l);
		fputs(header, tx_fp);
		for (i = 0; i < tx->l; i++) {
			if (i % LINELEN == 0)
				fputc('\n', tx_fp);
			fputc(tx->s[i], tx_fp);
		}
		fputc('\n', tx_fp);
	}
}

/**
 * Write one RNA-PET to the file
 *
 * rnapet: the file to store RNA-PETs;
 * tx: the transcript sequence.
 */
void spet(FILE *rnapet, const seq *tx, gene *gn) {
	char c, *pets = calloc(pet_len * 2 + 4, sizeof(char)), *header = malloc(
			BUFSIZE);
	int rand, i, j, shift, shift_2;

	// Simulate the RNA-PETs for a few copies
	rand = (int) (rand_f() * 4) + 1; // [1, 4]
	//	printf("Rounds = %d\n", rand);
	n_pets += rand * 2;
	for (j = 0; j <= rand; j++) {
		sprintf(header, ">%d\t %s\n", pet_id++, gn->id);
		fputs(header, rnapet);
		// RNA-PET introduces some shift
		shift = (int) (rand_f() * 4);
		shift_2 = (int) (rand_f() * 4);
		for (i = 0; i < pet_len; i++) {
			// From the beginning of the transcript
			c = tx->s[i + shift];
			if (rand_f() <= err) {
				c = mutate(c);
			}
			pets[i] = tolower(c);

			// From the end of transcript
			c = tx->s[tx->l - pet_len + i - shift_2];
			// Introduce some random error
			if (rand_f() <= err) {
				c = mutate(c);
			}
			pets[pet_len + i + 2] = to_upper_lower(c);
		}
		pets[pet_len] = '\n';
		pets[pet_len + 1] = '\0';
		fputs(pets, rnapet);
		sprintf(header, ">%d\t %s\n", pet_id++, gn->id);
		fputs(header, rnapet);
		pets[pet_len * 2 + 2] = '\n';
		pets[pet_len * 2 + 3] = '\0';
		fputs(&pets[0] + pet_len + 2, rnapet);
	}
}

/**
 * Simulate RNA-seq libraries
 */
void sseq(FILE *rnaseq, FILE *rnaseq_2, const seq *tx, const gene *gn) {
	int i, j, read_c, valid_rc, index_s, index_e, span;
	char c, *seqs = calloc(read_len * 2 + 4, sizeof(char));
	char *header = malloc(BUFSIZE);
	// 50% to save the reverse complement
	int is_rev_cmp = (int) (rand_f() * 2);

	// Simulate the RNA-Seq library.
	// Calculate the # of reads given coverage, length and read length.
	read_c = tx->l * coverage / (read_len * 2) + 1;
	//	fprintf(stderr, "%d reads are being generated...\n", read_c);
	j = 0;
	valid_rc = 0;

	for (i = 0; i < read_c; i++) {
		index_s = (int) (rand_f() * tx->l);
		span = normal(mu, sigma);
		index_e = index_s + span - 1;
		is_rev_cmp = (int) (rand_f() * 2);

		// If the random index is invalid, try at most 10 times
		if ((index_e + read_len >= tx->l) || (index_s < 0) || (index_e < 0)
				|| ((span + read_len * 2) > tx->l)) {
			// fprintf(stderr, "span = %d, tx->l = %zd\n", span, tx->l);
			if (j++ > 50) {
				j = 0; // Reset the counter, trying to generate the next read
			} else {
				i--; // Try again.
			}
			continue;
		}
		n_seqs += 2;
		//		sprintf(header, ">%d.1\t|%s|l=%d|r=%d|s=%d\n", seq_id, tx->h,
		//				index_s, index_e, span);
		// Copy character by character, allowing some error rate
		for (j = 0; j <= read_len - 1; j++) {
			c = tx->s[index_s + j];
			if (rand_f() <= err) {
				c = mutate(c);
			}
			seqs[j] = to_upper_lower(c);

			c = tx->s[index_e + j];
			if (rand_f() <= err) {
				c = mutate(c);
			}
			seqs[read_len + j + 2] = to_upper_lower(c);
		}
		seqs[read_len] = '\n';
		seqs[read_len + 1] = '\0';
		seqs[read_len * 2 + 2] = '\n';
		seqs[read_len * 2 + 3] = '\0';

		if (!strand_specific && is_rev_cmp) {
			sprintf(header, ">%d\t%s\t%d\t-\n", seq_id, gn->id, tx->l - read_len - index_e);
			fputs(header, rnaseq);
			rev_cmp(&seqs[0] + read_len + 2, read_len);
			fputs(&seqs[0] + read_len + 2, rnaseq);
			if (!both_files)
				seq_id++;
			sprintf(header, ">%d\t%s\t%d\t-\n", seq_id++, gn->id, tx->l - read_len - index_s);
			fputs(header, rnaseq_2);
			rev_cmp(seqs, read_len);
			fputs(&seqs[0], rnaseq_2);
		} else {
			sprintf(header, ">%d\t%s\t%d\n", seq_id, gn->id, index_s);
			fputs(header, rnaseq);
			fputs(&seqs[0], rnaseq);
			if (!both_files)
				seq_id++;
			sprintf(header, ">%d\t%s\t%d\n", seq_id++, gn->id, index_e);
			fputs(header, rnaseq_2);
			rev_cmp(&seqs[0] + read_len + 2, read_len);
			fputs(&seqs[0] + read_len + 2, rnaseq_2);
		}

		valid_rc++;
		if (span < min_l)
			min_l = span;
		if (span > max_l)
			max_l = span;
		total_l += span;
	}
	// fprintf(stderr, "Transcript length: %zd, %d reads generated!\n", tx->l,
	//		valid_rc);
}

/**
 * Simulate the transcripts, RNA-PETs and RNA-SEQs.
 *
 * gene_fp: list of genes;
 * fa: sequence of original sequence.
 */
void simulate(FILE *gene_fp, seq *fa) {
	char buf[BUFSIZ];
	char *attr[18], lib[BUFSIZE] = "", lib_2[BUFSIZE] = "";
	seq *tx = (seq*) malloc(sizeof(seq));
	gene *gn;
	int counter = 0, i = 0, len = 0;
	int si_pet = 1, si_tx = 1; // Flags indicating whether to simulate transcripts and RNA-PETs
	FILE *tx_fp;
	FILE *rnaseq, *rnaseq_2 = 0, *rnapet;
	all_len = g_array_new(FALSE, FALSE, sizeof(int));

	// To check whether need to simulate the transcripts:
	// If there is a parameter -p, do it; otherwise, if default file exists, do it; otherwise, omit.
	if (tx_fn != NULL) {
		tx_fp = xopen(tx_fn, "w");
		si_tx = 1;
	} else {
		tx_fp = fopen("read/tx.fa", "r");
		if (tx_fp != NULL) {
			si_tx = 0;
		} else {
			tx_fp = xopen("read/tx.fa", "w");
			si_tx = 1;
		}
	}

	// To check whether need to simulate the RNA-PET:
	// If there is a parameter -p, do it; otherwise, if default file exists, do it; otherwise, omit.
	if (rnapet_fn != NULL) {
		rnapet = xopen(rnapet_fn, "w");
		si_pet = 1;
	} else {
		rnapet = fopen("read/rnapet.fa", "r");
		if (rnapet != NULL) {
			si_pet = 0;
		} else {
			rnapet = xopen("read/rnapet.fa", "w");
			si_pet = 1;
		}
	}

	// Have to generate RNA-SEQs.
	if (both_files) {
		sprintf(lib, "read/rnaseq_%d_1.fa", mu);
		rnaseq = xopen(lib, "w");
		sprintf(lib_2, "read/rnaseq_%d_2.fa", mu);
		rnaseq_2 = xopen(lib_2, "w");
	} else {
		sprintf(lib, "read/rnaseq_%d.fa", mu);
		rnaseq = xopen(lib, "w");
	}

	tx->h = (char*) malloc(BUFSIZE);

	while (fgets(buf, sizeof(buf), gene_fp)) {
		// Omit the first line
		if (counter == 0) {
			counter++;
			continue;
		}
		if (counter % 400 == 0) {
			fprintf(stderr, "[simulate] %d genes have been simulated...\n",
					counter);
		}

		i = 0;
		gn = (gene*) malloc(sizeof(gene));
		n_tx++;

		// All fields from the gene list
		attr[0] = strtok(buf, "\t");
		while (attr[i] != NULL) { //ensure a pointer was found
			//			printf("fields[%d] = %s\n", i, fields[i]);
			attr[++i] = strtok(NULL, "\t"); //continue to tokenize the string
		}

		gn->id = attr[1];
		gn->ref = attr[2];
		gn->strand = attr[3][0];
		gn->exon_c = atoi(attr[8]);
		gn->exon_s = attr[9];
		gn->exon_e = attr[10];

		// fprintf(stderr, "Gene %s: ", gn->id);

		// Simulate the transcripts
		stx(tx_fp, tx, fa, gn, si_tx);
		// Simulate the RNA-PETs
		if (si_pet)
			spet(rnapet, tx, gn);
		// Simulate the RNA-SEQs
		if (both_files)
			sseq(rnaseq, rnaseq_2, tx, gn);
		else
			sseq(rnaseq, rnaseq, tx, gn);
		if (counter++ >= GENE_N)
			break;
	}
	ave_l = total_l / (n_seqs / 2);
	g_array_sort(all_len, (GCompareFunc) cmp_lens);
	for (i = 0; i < all_len->len; i++) {
		len += g_array_index(all_len, int, i);
		if (len > total_tx_len * 0.5) {
			opt_n50 = g_array_index(all_len, int, i);
			break;
		}
	}

	fprintf(stderr, "Library generated: min = %d, max = %d, ave = %d \n",
			min_l, max_l, ave_l);
	fprintf(stderr, "[simulate] # of exons: %d \n", n_exons);
	fprintf(stderr, "[simulate] # of genes: %d \n", n_tx);
	fprintf(stderr, "[simulate] Optimal N50 value: %d \n", opt_n50);
	fprintf(stderr, "[simulate] # of RNA-seqs: %d \n", n_seqs);
	fprintf(stderr, "[simulate] # of RNA-pets: %d \n", n_pets);
	fclose(gene_fp);
	fclose(rnapet);
	fclose(rnaseq);
	if (both_files)
		fclose(rnaseq_2);
	fclose(tx_fp);
}

int peta_simu(int argc, char *argv[]) {
	seq *fasta;
	char *gene_list = 0;
	FILE *fp;
	char *ref_seq;
	clock_t t = clock();

	while ((c = getopt(argc, argv, "c:g:e:p:t:m:s:b:l:a:d:")) >= 0) {
		switch (c) {
		case 'c': // Coverage
			coverage = atoi(optarg);
			break;
		case 'g': // Gene list file
			gene_list = strdup(optarg);
			break;
		case 'e': // Error rate
			err = atof(optarg);
			break;
		case 'm': // Average value for normal distribution
			mu = atoi(optarg);
			break;
		case 'd': // Standard deviation value for normal distribution
			sigma = atoi(optarg);
			break;
		case 'p': // File name for RNA-pet
			rnapet_fn = optarg;
			break;
		case 't': // File name for transcripts
			tx_fn = optarg;
			break;
		case 'b':
			both_files = atoi(optarg);
			break;
		case 'l':
			read_len = atoi(optarg);
			break;
		case 'a':
			pet_len = atoi(optarg);
			break;
		case 's':
			strand_specific = atoi(optarg);
			break;
		default:
			return usage();
		}
	}

	{
		if (optind + 1 > argc)
			return usage();

		// The file of original sequence.
		ref_seq = argv[optind];

		if (gene_list == NULL) {
			gene_list = strdup(ref_seq);
			strcat(gene_list, ".txt");
		}

		fprintf(stderr, "Coverage = %d\n", coverage);
		fprintf(stderr, "Gene list = %s\n", gene_list);
		fprintf(stderr, "Error rate = %.2f\n", err);
		fprintf(stderr, "Average span = %d\n", mu);
		fprintf(stderr, "Standard deviation = %d\n", sigma);
		fprintf(stderr, "Reference genome = %s\n", ref_seq);
		fprintf(stderr, "Separate paired reads = %d\n", both_files);
	}

	{
		fprintf(stderr, "Loading sequence from file: %s\n", ref_seq);
		fasta = read_seq(ref_seq);
		fprintf(stderr, "Length of sequence loaded: %zd \n", fasta->l);
		//		fprintf(stderr, "FASTA: %s \n", fasta->s);
		fp = xopen(gene_list, "r");
		simulate(fp, fasta);
	}

	fprintf(stderr, "Simulation Done: %.2f sec\n", (float) (clock() - t)
			/ CLOCKS_PER_SEC);
	free(gene_list);
	return 0;
}
