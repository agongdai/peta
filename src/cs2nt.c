#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "bwtaln.h"
#include "stdaln.h"

/*
  Here is a delicate example. ref_nt=ATTAAC(RBRBG), read_cs=RBBOG. If we
  decode as ATTGAC(RBGOG), there are one color change and one nt change;
  if we decode as ATTAAC(RBRBG), there are two color changes.

  In DP, if color quality is smaller than COLOR_MM, we will use COLOR_MM
  as the penalty; otherwise, we will use color quality as the
  penalty. This means we always prefer two consistent color changes over
  a nt change, but if a color has high quality, we may prefer one nt
  change.

  In the above example, the penalties of the two types of decoding are
  q(B)+25 and q(B)+q(O), respectively. If q(O)>25, we prefer the first;
  otherwise the second. Note that no matter what we choose, the fourth
  base will get a low nt quality.
 */

#define COLOR_MM 19
#define NUCL_MM  25

static const int nst_ntnt2cs_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4 };

/*
  {A,C,G,T,N} -> {0,1,2,3,4}
  nt_ref[0..size]: nucleotide reference: 0/1/2/3/4
  cs_read[0..size-1]: color read+qual sequence: base<<6|qual; qual==63 for N
  nt_read[0..size]: nucleotide read sequence: 0/1/2/3 (returned)
  btarray[0..4*size]: backtrack array (working space)
 */
void cs2nt_DP(int size, const uint8_t *nt_ref, const uint8_t *cs_read, uint8_t *nt_read, uint8_t *btarray)
{
	int h[8], curr, last;
	int x, y, xmin, hmin, k;

	// h[0..3] and h[4..7] are the current and last best score array, depending on curr and last

	// recursion: initial value
	if (nt_ref[0] >= 4) memset(h, 0, sizeof(int) << 2);
	else {
		for (x = 0; x != 4; ++x) h[x] = NUCL_MM;
		h[nt_ref[0]] = 0;
	}
	// recursion: main loop
	curr = 1; last = 0;
	for (k = 1; k <= size; ++k) {
		for (x = 0; x != 4; ++x) {
			int min = 0x7fffffff, ymin = 0;
			for (y = 0; y != 4; ++y) {
				int s = h[last<<2|y];
				if ((cs_read[k-1]&0x3f) != 63 && cs_read[k-1]>>6 != nst_ntnt2cs_table[1<<x|1<<y])
					s += ((cs_read[k-1]&0x3f) < COLOR_MM)? COLOR_MM : (cs_read[k-1]&0x3f); // color mismatch
				if (nt_ref[k] < 4 && nt_ref[k] != x) s += NUCL_MM; // nt mismatch
				if (s < min) {
					min = s; ymin = y;
				}
			}
			h[curr<<2|x] = min; btarray[k<<2|x] = ymin;
		}
		last = curr; curr = 1 - curr; // swap
	}
	// back trace
	hmin = 0x7fffffff; xmin = 0;
	for (x = 0; x != 4; ++x) {
		if (h[last<<2|x] < hmin) {
			hmin = h[last<<2|x]; xmin = x;
		}
	}
	nt_read[size] = xmin;
	for (k = size - 1; k >= 0; --k)
		nt_read[k] = btarray[(k+1)<<2 | nt_read[k+1]];
}
/*
  nt_read[0..size]: nucleotide read sequence: 0/1/2/3
  cs_read[0..size-1]: color read+qual sequence: base<<6|qual; qual==63 for N
  tarray[0..size*2-1]: temporary array
 */
uint8_t *cs2nt_nt_qual(int size, const uint8_t *nt_read, const uint8_t *cs_read, uint8_t *tarray)
{
	int k, c1, c2;
	uint8_t *t2array = tarray + size;
	// get the color sequence of nt_read
	c1 = nt_read[0];
	for (k = 1; k <= size; ++k) {
		c2 = nt_read[k]; // in principle, there is no 'N' in nt_read[]; just in case
		tarray[k-1] = (c1 >= 4 || c2 >= 4)? 4 : nst_ntnt2cs_table[1<<c1 | 1<<c2];
		c1 = c2;
	}
	for (k = 1; k != size; ++k) {
		int q = 0;
		if (tarray[k-1] == cs_read[k-1]>>6 && tarray[k] == cs_read[k]>>6) {
			q = (int)(cs_read[k-1]&0x3f) + (int)(cs_read[k]&0x3f) + 10;
		} else if (tarray[k-1] == cs_read[k-1]>>6) {
			q = (int)(cs_read[k-1]&0x3f) - (int)(cs_read[k]&0x3f);
		} else if (tarray[k] == cs_read[k]>>6) {
			q = (int)(cs_read[k]&0x3f) - (int)(cs_read[k-1]&0x3f);
		} // else, q = 0
		if (q < 0) q = 0;
		if (q > 60) q = 60;
		t2array[k] = nt_read[k]<<6 | q;
		if ((cs_read[k-1]&0x3f) == 63 || (cs_read[k]&0x3f) == 63) t2array[k] = 0;
	}
	return t2array + 1; // of size-2
}

