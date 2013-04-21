// Starting point
Function kmer_ext:
	ordered_kmer = [list of kmers with frequency >1, ordered by frequency decreasingly]
	hash_table = {the kmer hash table built before: kmer_integer->frequency}
	all_templates = []
	foreach kmer in ordered_kmer:
		new_template = kmer									// Start from the kmer
		template_ext(new_template, hash_table, kmer, 0)		// Extend to the right
		template_ext(new_template, hash_table, kmer, 1)		// Extend to the left
		if new_template.length > 50:
			all_templates.add(new_template)
	return all_templates

// Extend a template until no next char is available
Function template_ext(template, hash_table, query, direction):
	while True:
		c = get_next_char_by_kmers(hash_table, query, direction)
		if c == -1:
			print "Not extended any more, stop"
			break
		else:
			mark_kmer_used(hash_table, query)				// For every kmer used during this extension, mark the frequency as 0.
			extend_template(template, c, direction)
			shift_query(query, direction)
			

// Get next char, which has max frequency
Function get_next_char_by_kmers(hash_table, query, direction):
	counters = []
	next_probable_kmer = 0									// Kmers are stored as int64
	for (i = 0; i < 4; i++):
		counters[i] = 0
		next_probable_kmer = shift_bit(query, direction)
		counters[i] = get_kmer_count(hash_table, next_probable_kmer)
		// Get its reverse complement as well
		next_probable_kmer = rev_com_kmer(next_probable_kmer)
		counters[i] += get_kmer_count(hash_table, next_probable_kmer)
	// Count max frequency
	for (i = 0; i < 4; i++):
		max = (counters[i] > max) ? counters[i] : max;
	// If no probable next char, return -1
	if (max == 0)
		return -1
	// Return next char, which has most frequency
	for (i = 0; i < 4; i++):
		if (max == counters[i]):
			return i

Function mark_kmer_used(hash_table, query):
	hash_table[query] = 0
	hash_table[rev_com_kmer(query)] = 0