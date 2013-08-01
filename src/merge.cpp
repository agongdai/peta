/*
 * merge.cpp
 *
 *  Created on: Jul 29, 2013
 *      Author: carl
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "utils.h"
#include "tpl.hpp"
#include "peseq.h"
#include "bwtaln.h"
#include "junction.hpp"
#include "merge.hpp"
#include "k_hash.h"

/**
 * Merge the right template to the left template
 * The right template is not destoryed (its attributes are not altered), it is simply marked as 'DEAD'.
 * The rigth template is destoryed at hash_to_array in ass.cpp.
 * If param rev_com is 1, will merge the reverse complement of right template to left.
 */
int merge_tpls(tpl *left, tpl *right, int ol, int rev_com) {
	bwa_seq_t *r = NULL;
	int i = 0, connected = 0;
	int new_locus = 0, l_len = 0;
	junction *jun = NULL;
	ubyte_t c = 0;
	if (!left->alive || !right->alive || left->r_tail || right->l_tail) {
		show_debug_msg(
				__func__,
				"[WARNING] Merging templates not alive or with wrong ends: [%d, %d] and [%d, %d] \n",
				left->id, left->len, right->id, right->len);
		return 0;
	}

	// If they have junctions, just ignore
	if (tpls_have_junction(left, right)) {
		show_debug_msg(__func__,
				"[WARNING] These two templates have some junction, cannot be merged. \n");
		show_debug_msg(__func__, "Left  template: [%d, %d] \n", left->id,
				left->len);
		show_debug_msg(__func__, "RIght template: [%d, %d] \n", right->id,
				right->len);
		return 0;
	}

	// If trying to reverse complement it, there cannot be any connection on right
	if (rev_com) {
		if ((right->m_juncs && right->m_juncs->len > 0) || (right->b_juncs
				&& right->b_juncs->len > 0)) {
			show_debug_msg(
					__func__,
					"[WARNING] The right template is connected, cannot be reverse complement merged. \n");
			show_debug_msg(__func__, "Left  template: [%d, %d] \n", left->id,
					left->len);
			show_debug_msg(__func__, "RIght template: [%d, %d] \n", right->id,
					right->len);
			return 0;
		}
	}

	show_debug_msg(
			__func__,
			"Merging templates [%d, %d] and [%d, %d] with reverse complement: %d \n",
			left->id, left->len, right->id, right->len, rev_com);

    //p_tpl(left);
    //p_tpl(right);

	l_len = left->len;
	// Copy the template sequence from right to left
	for (i = ol; i < right->len; i++) {
		c = rev_com ? right->ctg->rseq[i] : right->ctg->seq[i];
		ext_con(left->ctg, c, 0);
	}
	left->len = left->ctg->len;
	set_rev_com(left->ctg);

	// Copy the tail of right to left
	if (right->r_tail) {
		left->r_tail = new_seq(right->r_tail, right->r_tail->len, 0);
	}

	// Add the reads on right template to the left template
	for (i = 0; i < right->reads->len; i++) {
		r = (bwa_seq_t*) g_ptr_array_index(right->reads, i);
		if (r->status == USED && r->contig_id == right->id) {
			// Update the locus
			new_locus = rev_com ? (right->len - r->contig_locus - r->len
					+ l_len - ol) : (r->contig_locus + l_len - ol);
			// If need to reverse complement the right template,
			//	need to reverse the reverse complement flag on the reads
			if (rev_com) {
				r->rev_com = r->rev_com ? 0 : 1;
			}
			add2tpl(left, r, new_locus);
		}
	}
	while (right->reads->len > 0)
		g_ptr_array_remove_index_fast(right->reads, 0);

    // The junctions on the right now go to the left template
	if (right->b_juncs) {
		for (i = 0; i < right->b_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(right->b_juncs, i);
			if (jun->branch_tpl == right) {
				jun->branch_tpl = left;
                // If right template connects to somewhere on the right, 
                //  now the left templates connects to it.
                if (!left->b_juncs)
                    left->b_juncs = g_ptr_array_sized_new(2);
                g_ptr_array_add(left->b_juncs, jun);
			}
		}
	}
	if (right->m_juncs) {
		for (i = 0; i < right->m_juncs->len; i++) {
			jun = (junction*) g_ptr_array_index(right->m_juncs, i);
			if (jun->main_tpl == right) {
				jun->main_tpl = left;
                // Update the junction locus
				jun->locus += l_len;
                if (!left->m_juncs)
                    left->m_juncs = g_ptr_array_sized_new(2);
                g_ptr_array_add(left->m_juncs, jun);
			}
		}
	}

    //p_tpl(left);
	right->alive = 0;
    return 1;
}
