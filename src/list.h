/* mystack.h -- Stack declaration and function prototypes:  */
#ifndef LIST_H_
#define LIST_H_

#include <glib.h>
#include "kmer.h"

typedef struct
{
	GPtrArray *starts;
	GPtrArray *kmers;
	int order; 			// 0: increasing; 1: reversed; -1: not sorted
} slist;

slist *new_slist();
void free_slist(slist *sl);
int slist_binary(mer *kmers, const uint32_t size, uint64_t s);
mer *slist_ins_pt(mer *kmers, uint32_t *size, uint32_t *full_space, mer *new_mer);
gint cmp_kmer_by_seq(gpointer a, gpointer b);
gint cmp_kmer_by_count(gpointer a, gpointer b);

#endif
