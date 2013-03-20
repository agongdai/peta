/* mystack.h -- Stack declaration and function prototypes:  */
#ifndef LIST_H_
#define LIST_H_

#include <glib.h>
#include "bwase.h"

typedef struct
{
	GPtrArray *values;
	GPtrArray *kmers;
	int order; 			// 0: increasing; 1: reversed; -1: not sorted
} slist;

slist *new_slist();
void free_slist(slist *sl);
int slist_binary(slist *sl, bwa_seq_t *query);
int slist_ins_pt(slist *sl, bwa_seq_t *new_pt);
gint cmp_kmer_by_seq(gpointer a, gpointer b);
gint cmp_kmer_by_freq(gpointer a, gpointer b);

#endif
