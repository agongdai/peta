#include <vector>
#include <unordered_map>
#include <cassert>
#include <unistd.h>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <stdint.h>
#include <glib.h>
#include <time.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <string.h>
#include <stdio.h>
#include "kmers.h"
#include "utils.h"
#include "pehash.h"
#include "list.h"
#include "peseq.h"
#include "rnaseq.h"
#include "bwtaln.h"
#include "edge.h"
#include "edgelist.h"
#include "pool.h"
#include "mate.h"
#include "merge.h"
#include "readrm.h"

using namespace std;

int kmer_ctg_id = 0;
int ins_size = 0;
int sd_ins_size = 0;
int kmer_n_threads = 0;
int kmer_len = 0;
char *kmer_out = NULL;
GMutex *kmer_id_mutex;
struct timespec kmer_start_time, kmer_finish_time;

mer_meta *new_mer_meta() {
	mer_meta *m = (mer_meta*) malloc(sizeof(mer_meta));
	m->n_seqs = 0;
	m->n_kmers = 0;
	m->k = 0;
	return m;
}

mer *new_mer() {
	mer *m = (mer*) malloc(sizeof(struct mer));
	m->count = 0;
	m->s = 0;
	m->status = FRESH;
	return m;
}

void free_kmer_list(GPtrArray *kmer_list) {
	uint64_t i = 0;
	mer *m = NULL;
	for (i = 0; i < kmer_list->len; i++) {
		m = (mer*) g_ptr_array_index(kmer_list, i);
		free(m);
	}
	g_ptr_array_free(kmer_list, TRUE);
}

uint64_t get_kmer_int(const ubyte_t *seq, const int start,
		const int interleaving, const int len) {
	uint64_t key = 0;
	int i = 0;
	for (i = 0; i < len * interleaving; i += interleaving) {
		if (seq[i + start] < 0 || seq[i + start] > 4) {
			return 0;
		}
		key <<= 2;
		// printf("%d \n", seq[i + start]);
		// Keep the last two digits and perform 'or' oepration
		key = key | (3 & seq[i + start]);
		// printf("%d: key = %" ID64 "\n", i, key);
	}
	return key;
}

void build_kmers_hash(const char *fa_fn, hash_opt *opt) {
	clock_t t = clock();
	uint32_t n_reads = 0, i = 0, j = 0;
	uint64_t mer_v = 0, *kmer_freq = NULL;
	mer_counter counter;
	mer_hash hash;
	bwa_seq_t *reads = NULL, *r = NULL;
	FILE *hash_fp;
	char *hash_fn = (char*) malloc(BUFSIZE);

	show_msg(__func__, "Hashing library %s...\n", fa_fn);
	reads = load_reads(fa_fn, &n_reads);

	// Count frequencies of every kmer
	show_msg(__func__, "Counting kmer frequencies ...\n");
	for (i = 0; i < n_reads; i++) {
		r = &reads[i];
		for (j = 0; j <= r->len - opt->k; j++) {
			mer_v = get_kmer_int(r->seq, j, 1, opt->k);
			if (counter[mer_v]) {
				counter[mer_v]++;
			} else {
				counter[mer_v] = 1;
				opt->n_k_mers++;
			}
			opt->n_pos++;
		}
	}

	// Set the value of hash map: first value stores how many reads having this kmer
	show_msg(__func__, "Allocating space for %" ID64 " kmers, %" ID64 " distinct kmers ...\n", opt->n_pos, opt->n_k_mers);
	for (mer_counter::iterator it = counter.begin(); it != counter.end(); ++it) {
		kmer_freq = (uint64_t*) calloc(it->second + 1, sizeof(uint64_t));
		kmer_freq[0] = 0;
		hash[it->first] = kmer_freq;
	}
	show_msg(__func__, "Setting hash map ...\n");
	for (i = 0; i < n_reads; i++) {
		r = &reads[i];
		for (j = 0; j <= r->len - opt->k; j++) {
			mer_v = get_kmer_int(r->seq, j, 1, opt->k);
			kmer_freq = hash[mer_v];
			kmer_freq[++kmer_freq[0]] = get_hash_value(atoi(r->name), j);
		}
	}
	show_msg(__func__, "Saving hash map of %d reads: %.2f sec ... \n", n_reads,
			(float) (clock() - t) / CLOCKS_PER_SEC);
	sprintf(hash_fn, "%s.map", fa_fn);
	hash_fp = xopen(hash_fn, "w");
	fwrite(opt, sizeof(hash_opt), 1, hash_fp);
	for (mer_hash::iterator m = hash.begin(); m != hash.end(); ++m) {
		mer_v = m->first;
		fwrite(&mer_v, sizeof(uint64_t), 1, hash_fp);
		kmer_freq = m->second;
		fwrite(kmer_freq, sizeof(uint64_t), kmer_freq[0] + 1, hash_fp);
	}
	show_msg(__func__, "%d reads hashed: %.2f sec ... \n", n_reads,
			(float) (clock() - t) / CLOCKS_PER_SEC);
}

hash_map *load_hash_map(const char *fa_fn, mer_hash& kmers) {
	clock_t t = clock();
	FILE *map_fp = NULL;
	hash_opt *opt = NULL;
	uint64_t i = 0, *count = NULL, *freq = NULL, *kmer_int = NULL;
	uint32_t n_seqs = 0;
	int rs = 0;
	bwa_seq_t *seqs = NULL;
	char *map_fn = (char*) malloc(BUFSIZE);
	hash_map *hm = (hash_map*) malloc(sizeof(hash_map));

	seqs = load_reads(fa_fn, &n_seqs);
	hm->seqs = seqs;
	hm->n_seqs = n_seqs;

	sprintf(map_fn, "%s.map", fa_fn);
	show_msg(__func__, "Loading hash map from %s ...\n", map_fn);
	map_fp = xopen(map_fn, "rb");
	opt = (hash_opt*) malloc(sizeof(hash_opt));
	if (!fread(opt, sizeof(hash_opt), 1, map_fp)) {
		err_fatal(__func__, "Unable to read from the map file %s! \n", map_fn);
	}
	for (i = 0; i < opt->n_k_mers; i++) {
		kmer_int = (uint64_t*) malloc(sizeof(uint64_t));
		rs = fread(kmer_int, sizeof(uint64_t), 1, map_fp);
		count = (uint64_t*) malloc(sizeof(uint64_t));
		rs = fread(count, sizeof(uint64_t), 1, map_fp);
		freq = (uint64_t*) calloc(*count + 1, sizeof(uint64_t));
		rs = fread(freq + 1, sizeof(uint64_t), *count, map_fp);
		freq[0] = *count;
		kmers[*kmer_int] = freq;
	}

	fclose(map_fp);
	free(map_fn);
	hm->o = opt;
	hm->hash = &kmers;
	show_msg(__func__, "%" ID64 " distinct kmers loaded\n", opt->n_k_mers);
	show_msg(__func__, "%" ID64 " kmer occurrences loaded\n", opt->n_pos);
	show_msg(__func__, "Hash map loaded: %.2f sec ... \n",
			(float) (clock() - t) / CLOCKS_PER_SEC);
	return hm;
}

void test_kmer_hash(const char *fa_fn) {
	mer_hash map, *map_copy = NULL;
	uint64_t i = 0, *freq = NULL;
	uint64_t j = 0, seq_id = 0;
	int pos = 0;
	hash_map *hm = load_hash_map(fa_fn, map);
	map_copy = hm->hash;
	for (i = 0; i < 5; i++) {
		freq = (*map_copy)[i];
		for (j = i; j < freq[0]; j++) {
			read_hash_value(&seq_id, &pos, freq[j]);
show_msg		(__func__, "Read %" ID64 " contains kmer %" ID64 " at %d \n", seq_id, i, pos);
	}
}
}

void build_kmers(const char *fa_fn, const char *out_fn, const int k) {
	bwa_seq_t *reads = NULL, *r = NULL;
	uint32_t n_reads = 0, i = 0;
	int j = 0;
	uint64_t mer_v = 0, n_kmers = 0;
	mer *m = NULL;
	mer_meta *meta = NULL;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);

	std::unordered_map<uint64_t, mer*> hash;
	GPtrArray* kmer_list = g_ptr_array_sized_new(BUFSIZ);
	FILE *out = NULL;
	show_msg(__func__, "Building the kmers hashtable of %s... \n", fa_fn);
	reads = load_reads(fa_fn, &n_reads);
	for (i = 0; i < n_reads; i++) {
		r = &reads[i];
		//p_query(__func__, r);
		if (i % 1000000 == 0)
			show_debug_msg(__func__,
					"%d sequences kmer counted: %d kmers...\n", i, n_kmers);
		for (j = 0; j < r->len - k; j++) {
			mer_v = get_kmer_int(r->seq, j, 1, k);
			//show_debug_msg(__func__, "Kmer id: %" ID64 " \n", *mer_v);
			m = hash[mer_v];
			if (m) {
				m = (mer*) m;
				m->count++;
				//show_debug_msg(__func__, "Kmer exists: %" ID64 ", %d \n", exist_m->s, exist_m->count);
			} else {
				m = new_mer();
				m->s = mer_v;
				m->count++;
				n_kmers++;
				//show_debug_msg(__func__, "Kmer added: %" ID64 ", %d \n", m->s, m->count);
				g_ptr_array_add(kmer_list, m);
				hash[mer_v] = m;
			}
		}
	}
	meta = new_mer_meta();
	meta->k = k;
	meta->n_kmers = n_kmers;
	meta->n_seqs = n_reads;
	show_debug_msg(__func__, "%d distinct %d-mers counted \n", n_kmers, k);
	g_ptr_array_sort(kmer_list, (GCompareFunc) cmp_kmer_by_count);
	out = xopen(out_fn, "w");
	fwrite(meta, sizeof(mer_meta), 1, out);
	for (i = 0; i < kmer_list->len; i++) {
		m = (mer*) g_ptr_array_index(kmer_list, i);
		fwrite(m, sizeof(mer), 1, out);
	}
	free_kmer_list(kmer_list);
	bwa_free_read_seq(n_reads, reads);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	fclose(out);
}

void load_kmers(const char *kmer_file, mer_map& kmers, GPtrArray *kmer_list,
		mer_meta *meta) {
	FILE *kmers_f = xopen(kmer_file, "r");
	uint64_t i = 0, reading = 0;
	mer *m = NULL;

	show_msg(__func__, "Loading kmers...\n");
	reading = fread(meta, sizeof(mer_meta), 1, kmers_f);
	show_msg(__func__, "%d %d-mers to load... \n", meta->n_kmers, meta->k);
	for (i = 0; i < meta->n_kmers; i++) {
		if (i && i % 10000000 == 0)
			show_msg(__func__, "%d kmers loaded...\n", i);
		m = new_mer();
		reading = fread(m, sizeof(struct mer), 1, kmers_f);
		kmers[m->s] = m;
		g_ptr_array_add(kmer_list, m);
	}
	fclose(kmers_f);
}

bwa_seq_t *get_kmer_seq(mer *kmer, const int k) {
	ubyte_t *seq = NULL;
	bwa_seq_t *read = NULL;
	int i = 0;
	uint64_t copy = kmer->s, index = 0;
	read = blank_seq();
	free(read->seq);
	read->seq = (ubyte_t*) calloc(k + 1, sizeof(ubyte_t));
	seq = read->seq;
	for (i = 0; i < k; i++) {
		index = copy;
		index &= 3; // Keep only the last two bits
		seq[k - 1 - i] = index;
		copy >>= 2;
	}
	read->len = k;
	read->name = (char*) malloc(32);
	sprintf(read->name, "%d", kmer->count);
	set_rev_com(read);
	return read;
}

/**
 * Check kmer frequencies of all possible next chars, get the highest one
 * Return: ACGTN -> 01234;
 *         None -> -1;
 */
int next_char_by_kmers(mer_map *kmers, mer_meta *meta, bwa_seq_t *query,
		const int ori) {
	bwa_seq_t *next_q = NULL;
	int counters[5], i = 0, max = 0;
	uint64_t kmer_int = 0;
	mer *m = NULL;
	for (i = 0; i < 5; i++) {
		next_q = new_seq(query, query->len, 0);
		ext_que(next_q, i, ori);
		counters[i] = 0;
		// Check the forward kmer
		//p_query("SEQ", next_q);
		kmer_int = get_kmer_int(next_q->seq, 0, 1, meta->k);
		m = (*kmers)[kmer_int];
		if (m) {
			//show_debug_msg(__func__, "Count: %d\n", m->count);
			if (m->status != USED)
				counters[i] = m->count;
		}
		// Check the reverse kmer
		kmer_int = get_kmer_int(next_q->rseq, 0, 1, meta->k);
		m = (*kmers)[kmer_int];
		//p_query("REV", next_q);
		if (m) {
			//show_debug_msg(__func__, "Count: %d\n", m->count);
			if (m->status != USED)
				counters[i] += m->count;
		}
		bwa_free_read_seq(1, next_q);
	}
	// Get max occurrence of the chars
	for (i = 0; i < 5; i++) {
		max = (counters[i] > max) ? counters[i] : max;
	}
	//show_debug_msg(__func__, "Max: %d \n", max);
	if (max == 0)
		return -1;
	for (i = 0; i < 5; i++) {
		if (max == counters[i])
			return i;
	}
	return -1;
}

/**
 * Mark a kmer as used.
 * The parameter 'kmer' is a 'copy' only, need to mark the kmer on the list
 */
void mark_kmer_used(mer_map *kmers, mer_meta *meta, bwa_seq_t *query) {
	uint64_t kmer_int = 0;
	mer *m = NULL;
	kmer_int = get_kmer_int(query->seq, 0, 1, meta->k);
	m = (*kmers)[kmer_int];
	if (m) {
		m->status = USED;
	}
	kmer_int = get_kmer_int(query->rseq, 0, 1, meta->k);
	m = (*kmers)[kmer_int];
	if (m) {
		m->status = USED;
	}
}

int ext_eg_query_n(mer_map *kmers, mer_meta *meta, bwa_seq_t *query, edge *eg,
		bwa_seq_t *read, const int shift, const int ori) {
	int i = 0, ext_len = 0;
	ubyte_t *s = NULL, next_c = 0;
	if (shift > query->len)
		return 0;
	s = read->seq;
	for (i = shift; i < read->len; i++) {
		next_c = s[i];
		ext_que(query, next_c, ori);
		ext_con(eg->contig, next_c, 0);
		//p_ctg_seq(__func__, eg->contig);
		mark_kmer_used(kmers, meta, query);
	}
	ext_len = eg->contig->len - eg->len;
	eg->len = eg->contig->len;
	return ext_len;
}

int ext_by_mates(mer_map *kmers, mer_meta *meta, edgearray *mates,
		bwa_seq_t *query, edge *eg, const int ori) {
	int ol = 0, ext_len = 0;
	uint32_t i = 0;
	bwa_seq_t *m = NULL, *tmp = NULL, *tpl = NULL;
	tpl = new_seq(eg->contig, query->len, eg->len - query->len);
	if (ori) {
		seq_reverse(tpl->len, tpl->seq, 0);
	}
	for (i = 0; i < mates->len; i++) {
		m = (bwa_seq_t*) g_ptr_array_index(mates, i);
		tmp = m;
		if (m->rev_com)
			tmp = new_mem_rev_seq(m, m->len, 0);
		ol = find_ol_within_k(tmp, tpl, 1, 0, tpl->len, ori);
		//p_ctg_seq(__func__, tmp);
		//p_ctg_seq(__func__, template);
		//show_debug_msg(__func__, "OL: %d \n", ol);
		if (ol >= 8) {
			//p_query("Mate", m);
			ext_len = ext_eg_query_n(kmers, meta, query, eg, tmp, ol, ori);
		}

		if (m->rev_com)
			bwa_free_read_seq(1, tmp);
		if (ext_len > 0)
			break;
	}
	bwa_free_read_seq(1, tpl);
	return ext_len;
}

void kmer_ext_edge(edge *eg, bwa_seq_t *query, mer_map *kmers, mer_meta *meta,
		hash_table *ht, const int ori) {
	int next_c = 0, pre_round_len = eg->len;
	pool *mate_pool = NULL;
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	p_query(__func__, query);
	while (1) {
		//p_query(__func__, query);
		//p_ctg_seq("Contig", eg->contig);
		next_c = next_char_by_kmers(kmers, meta, query, ori);
		//show_debug_msg(__func__, "Next char: %d \n", next_c);
		// If cannot find next char
		if (next_c == -1) {
			show_debug_msg(__func__, "[%d, %d] Realigning reads to edge \n",
					eg->id, eg->len);
			//p_ctg_seq("Contig", eg->contig);
			realign_extended(ht, eg, pre_round_len, 0, ori);
			//p_readarray(eg->reads, 1);
			//show_debug_msg(__func__, "Getting mate pool \n");
			mate_pool = get_mate_pool_from_edge(eg, ht, ori, ins_size,
					sd_ins_size);
			//p_readarray(mate_pool->reads, 1);
			show_debug_msg(__func__, "[%d, %d] Extending by mates \n", eg->id,
					eg->len);
			if (ext_by_mates(kmers, meta, mate_pool->reads, query, eg, ori)) {
				pre_round_len = eg->len;
				continue;
			} else
				break;
		}
		ext_con(eg->contig, next_c, 0);
		eg->len = eg->contig->len;
		ext_que(query, next_c, ori);
		mark_kmer_used(kmers, meta, query);
		if (eg->len % 100 == 0)
			show_debug_msg(__func__, "Edge [%d, %d] \n", eg->id, eg->len);
	}
	if (ori)
		seq_reverse(eg->len, eg->contig->seq, 0);
	bwa_free_read_seq(1, query);
}

typedef struct {
	mer_map *kmers;
	mer_meta *meta;
	hash_table *ht;
	GPtrArray *edge_list;
} kmer_thread_meta;

/**
 * Extend a kmer
 */
void *kmer_ext_thread(gpointer data, gpointer thread_params) {
	edge *eg = NULL;
	mer *m = NULL;
	bwa_seq_t *query = NULL, *kmer = NULL;
	int round_1_len = 0, round_2_len = 0;
	kmer_thread_meta *params = (kmer_thread_meta*) thread_params;
	mer_map *kmers = (params->kmers);
	m = (mer*) data;
	kmer = get_kmer_seq(m, params->meta->k);

	if (has_n(kmer, 4) || is_biased_q(kmer)) {
		bwa_free_read_seq(1, kmer);
		return NULL;
	}show_debug_msg(__func__,
			"========= Kmer: %" ID64 " ========= \n", get_kmer_int(kmer->seq, 0, 1, kmer->len));
	eg = new_eg();
	g_mutex_lock(kmer_id_mutex);
	eg->id = kmer_ctg_id++;
	g_mutex_unlock(kmer_id_mutex);
	// Get a copy of the kmer
	query = new_seq(kmer, kmer->len, 0);
	eg->contig = new_seq(query, query->len, 0);
	eg->len = eg->contig->len;
	//	p_query(__func__, kmer);
	//	p_query(__func__, query);
	kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 0);
	round_1_len = eg->len;
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	query = new_seq(eg->contig, kmer->len, 0);
	kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 1);
	//p_ctg_seq("Contig", eg->contig);
	show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);

	if (eg->len - round_1_len > 2) {
		query = new_seq(eg->contig, kmer->len, 0);
		round_2_len = eg->len;
		kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 0);
		show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		if (eg->len - round_2_len > 2) {
			query = new_seq(eg->contig, kmer->len, 0);
			kmer_ext_edge(eg, query, kmers, params->meta, params->ht, 1);
			show_debug_msg(__func__, "Edge %d length: %d \n", eg->id, eg->len);
		}
	}

	eg->start_kmer_int = get_kmer_int(kmer->seq, 0, 1, kmer->len);
	if (eg->len < 100000)
		destroy_eg(eg);
	else {
		g_mutex_lock(kmer_id_mutex);
		g_ptr_array_add(params->edge_list, eg);
		g_mutex_unlock(kmer_id_mutex);
	}
	bwa_free_read_seq(1, kmer);

	return NULL;
}

void kmer_threads(kmer_thread_meta *params, GPtrArray *kmer_list) {
	GThreadPool *thread_pool = NULL;
	mer *m = NULL;
	uint64_t i = 0;
	thread_pool = g_thread_pool_new((GFunc) kmer_ext_thread, params, 1, TRUE,
			NULL);

	show_msg(__func__, "Extending by kmers...\n");
	for (i = 0; i < kmer_list->len; i++) {
		m = (mer*) g_ptr_array_index(kmer_list, i);
		if (m->status != USED && m->count > 1) {
			g_thread_pool_push(thread_pool, (gpointer) m, NULL);
		}
	}
	g_thread_pool_free(thread_pool, 0, 1);
}

void test_kmer_ext(kmer_thread_meta *params, GPtrArray *kmer_list) {
	FILE *kmer_contigs = NULL;
	bwa_seq_t *kmer_seq = NULL;
	mer_map *map = params->kmers;

	mer *m = (mer*) (*map)[53073660363574];
	kmer_seq = get_kmer_seq(m, params->meta->k);
	kmer_ext_thread((gpointer) kmer_seq, (gpointer) params);

	m = (mer*) (*map)[857693345222861];
	kmer_seq = get_kmer_seq(m, params->meta->k);
	kmer_ext_thread((gpointer) kmer_seq, (gpointer) params);

	kmer_contigs = xopen("../SRR097897_out/single.fa", "w");
	save_edges(params->edge_list, kmer_contigs, 0, 0, 0);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
}

void ext_by_kmers_core(char *lib_file, const char *solid_file,
		const char *kmer_file, const int insert_size, const int sd_insert_size,
		const int n_threads) {
	hash_table *ht = NULL;
	GPtrArray *all_edges = NULL, *kmer_list = NULL;
	mer_meta *meta = NULL;
	FILE *kmer_contigs = NULL;
	mer_map kmers;
	kmer_thread_meta *params = (kmer_thread_meta*) calloc(1,
			sizeof(kmer_thread_meta));

	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);

	show_msg(__func__, "Library: %s \n", lib_file);
	show_msg(__func__, "Solid Reads: %s \n", solid_file);
	ins_size = insert_size;
	sd_ins_size = sd_insert_size;

	kmer_list = g_ptr_array_sized_new(BUFSIZ);
	meta = (mer_meta*) malloc(sizeof(mer_meta));
	load_kmers(kmer_file, kmers, kmer_list, meta);
	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done loading kmers: %.2f sec\n",
			(float) (kmer_finish_time.tv_sec - kmer_start_time.tv_sec));

	ht = pe_load_hash(lib_file);
	all_edges = g_ptr_array_sized_new(BUFSIZ);
	params->ht = ht;
	params->kmers = &kmers;
	params->meta = meta;
	params->edge_list = all_edges;

	//test_kmer_ext(params, kmer_list);
	//exit(1);
	kmer_threads(params, kmer_list);

	kmer_contigs = xopen("../SRR097897_out/kmer_contigs.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
	g_ptr_array_free(kmer_list, TRUE);

	reset_edge_ids(all_edges);
	merge_ol_edges(all_edges, insert_size, sd_insert_size, ht, n_threads);
	destroy_ht(ht);
	free(params);

	kmer_contigs = xopen("../SRR097897_out/merged.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 100);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
	kmer_contigs = xopen("../SRR097897_out/merged_all.fa", "w");
	save_edges(all_edges, kmer_contigs, 0, 0, 0);
	fflush(kmer_contigs);
	fclose(kmer_contigs);
}

int pe_kmer(int argc, char *argv[]) {
	int c = 0;
	clock_gettime(CLOCK_MONOTONIC, &kmer_start_time);
	while ((c = getopt(argc, argv, "k:m:s:o:t:")) >= 0) {
		switch (c) {
		case 'k':
			kmer_len = atoi(optarg);
			break;
		case 'o':
			kmer_out = optarg;
			break;
		case 'm':
			ins_size = atoi(optarg);
			break;
		case 's':
			sd_ins_size = atoi(optarg);
			break;
		case 't':
			kmer_n_threads = atoi(optarg);
			break;
		}
	}
	if (optind + 2 > argc) {
		show_msg(__func__, "Parameters error! \n");
		return 1;
	}
	if (!g_thread_supported())
		g_thread_init(NULL);
	kmer_id_mutex = g_mutex_new();

	hash_opt *opt = (hash_opt*) malloc(sizeof(hash_opt));
	opt->k = kmer_len;
	test_kmer_hash(
			"/home/carl/Projects/peta/rnaseq/Spombe/SRR097897/SRR097897_corrected.fa");
	//build_kmers_hash(argv[optind], opt);
	//ext_by_kmers_core(argv[optind], argv[optind + 1], argv[optind + 2],
	//		ins_size, sd_ins_size, kmer_n_threads);

	clock_gettime(CLOCK_MONOTONIC, &kmer_finish_time);
	show_msg(__func__, "Done: %.2f sec\n", (float) (kmer_finish_time.tv_sec
			- kmer_start_time.tv_sec));
	return 0;
}
