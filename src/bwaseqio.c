#include <zlib.h>
#include <ctype.h>
#include "bwtaln.h"
#include "utils.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];

struct __bwa_seqio_t {
	// for BAM input
	int is_bam, which; // 1st bit: read1, 2nd bit: read2, 3rd: SE
	// for fastq input
	kseq_t *ks;
};

bwa_seqio_t *bwa_seq_open(const char *fn) {
	gzFile fp;
	bwa_seqio_t *bs;
	bs = (bwa_seqio_t*) calloc(1, sizeof(bwa_seqio_t));
	fp = xzopen(fn, "r");
	bs->ks = kseq_init(fp);
	return bs;
}

void bwa_seq_close(bwa_seqio_t *bs) {
	if (bs == 0)
		return;
	gzclose(bs->ks->f->f);
	kseq_destroy(bs->ks);
	free(bs);
}

void seq_reverse(int len, ubyte_t *seq, int is_comp) {
	int i;
	if (is_comp) {
		for (i = 0; i < len >> 1; ++i) {
			char tmp = seq[len - 1 - i];
			if (tmp < 4)
				tmp = 3 - tmp;
			seq[len - 1 - i] = (seq[i] >= 4) ? seq[i] : 3 - seq[i];
			seq[i] = tmp;
		}
		if (len & 1)
			seq[i] = (seq[i] >= 4) ? seq[i] : 3 - seq[i];
	} else {
		for (i = 0; i < len >> 1; ++i) {
			char tmp = seq[len - 1 - i];
			seq[len - 1 - i] = seq[i];
			seq[i] = tmp;
		}
	}
}

#define BARCODE_LOW_QUAL 13

bwa_seq_t *bwa_read_seq(bwa_seqio_t *bs, int n_needed, int *n, int mode,
		int trim_qual) {
	bwa_seq_t *seqs, *p;
	kseq_t *seq = bs->ks;
	int n_seqs, l, i = 0, is_comp = mode & BWA_MODE_COMPREAD, is_64 = mode
			& BWA_MODE_IL13, l_bc = mode >> 24;
	long n_trimmed = 0, n_tot = 0;

	if (l_bc > 15) {
		fprintf(stderr, "[%s] the maximum barcode length is 15.\n", __func__);
		return 0;
	}
	n_seqs = 0;
	seqs = (bwa_seq_t*) calloc(n_needed, sizeof(bwa_seq_t));
	while ((l = kseq_read(seq)) >= 0) {
		if (is_64 && seq->qual.l)
			for (i = 0; i < seq->qual.l; ++i)
				seq->qual.s[i] -= 31;
		if (seq->seq.l <= l_bc)
			continue; // sequence length equals or smaller than the barcode length
		p = &seqs[n_seqs++];
		if (l_bc) { // then trim barcode
			for (; i < seq->seq.l; ++i)
				seq->seq.s[i - l_bc] = seq->seq.s[i];
			seq->seq.l -= l_bc;
			seq->seq.s[seq->seq.l] = 0;
			if (seq->qual.l) {
				for (i = l_bc; i < seq->qual.l; ++i)
					seq->qual.s[i - l_bc] = seq->qual.s[i];
				seq->qual.l -= l_bc;
				seq->qual.s[seq->qual.l] = 0;
			}
			l = seq->seq.l;
		} //else
		p->full_len = p->len = l;
		n_tot += p->full_len;
		p->seq = (ubyte_t*) calloc(p->len, 1);

		p->contig_id = -1; // -1 means unused
		p->contig_locus = 0;
		p->rev_com = 0;
		p->pos = -1;
		p->cursor = -1;

		for (i = 0; i != p->full_len; ++i) {
			//			fprintf(stderr, "%c", seq->seq.s[i]);
			p->seq[i] = nst_nt4_table[(int) seq->seq.s[i]];
			//			fprintf(stderr, "%d", (int) seqs->seq[i]);
		}

		p->rseq = (ubyte_t*) calloc(p->full_len, 1);
		memcpy(p->rseq, p->seq, p->len);
		seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = strdup((const char*) seq->name.s);
		{ // trim /[12]$
			int t = strlen(p->name);
			if (t > 2 && p->name[t - 2] == '/' && (p->name[t - 1] == '1'
					|| p->name[t - 1] == '2'))
				p->name[t - 2] = '\0';
		}
		if (n_seqs == n_needed)
			break;
	}
	*n = n_seqs;
	if (n_seqs && trim_qual >= 1)
		fprintf(stderr, "[bwa_read_seq] %.1f%% bases are trimmed.\n", 100.0f
				* n_trimmed / n_tot);
	if (n_seqs == 0) {
		free(seqs);
		return 0;
	}
	return seqs;
}

void bwa_free_read_seq(int n_seqs, bwa_seq_t *seqs) {
	int i;
	if (seqs && n_seqs > 0) {
		for (i = 0; i != n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			free(p->name);
			free(p->seq);
			free(p->rseq);
		}
		free(seqs);
	}
}
