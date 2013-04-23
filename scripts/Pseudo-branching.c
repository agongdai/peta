// Global values
ALL_TEMPLATES = []			// Store all assembled templates
AMBIGUOUS_BASES = []		// Element is in format of (template, locus), just for probable error correction
BRANCHING_EVENTS = []       // Element is in format of (main_template, branch, locus, direction), for post-processing

// Starting point
Function kmer_ext:
	ordered_kmer = [list of kmers with frequency >1, ordered by frequency decreasingly]
	hash_table = {the kmer hash table built before: kmer_integer->[frequency, read_1, read_2...] }
	foreach kmer in ordered_kmer:
		new_template = kmer											// Start from the kmer
		template_ext(new_template, hash_table, kmer, 0)				// Extend to the right
		template_ext(new_template, hash_table, kmer, 1)				// Extend to the left
		update_branching_events(BRANCHING_EVENTS, new_template)		// The locus info needs to update
		update_ambiguous_bases(AMBIGUOUS_BASES, new_template)		// The locus info needs to update
		ALL_TEMPLATES.add(template)

/**
	A recursive function to extend templates
	Main process (together with Function kmer_ext):
		1. Start from some kmer;
		2. Extend the template to the right, based on frequencies of next probable kmers;
		3. If next kmer is ambiguous, and there are junction reads supporting branching:
			3.1 extend the main template (repeat step 2)
			3.2 extend the branch template (repeat step 2), to one direction only:
			3.3 record the branching event to BRANCHING_EVENTS
			3.4 add the branch template to ALL_TEMPLATES
		4. Extend the template to the left. If there is branching, repeat step 3 similarly
		5. Add current template to ALL_TEMPLATES
**/
Function template_ext(template, hash_table, query, direction):
	counters = []
	while True:
		counters = count_next_kmers(hash_table, query, direction)
		max_c = get_max_index(counters)
		if max_c == -1:
			print "Not extended any more, stop"
			break
		second_c = get_second_freq_char(counters, max_c)
		// If the second largest kmer-frequency is at least 1, check branching
		if (second_c != -1):
			// If it is just noise, likely to come to the same sequence quickly
			// So here try to extend extra 5 bases quickly, for both branches
			main_query = shift_query(query, max_c, direction)
			main_branch = try_short_tpl_ext(hash_table, main_query, direction, 5)
			second_query = shift_query(query, second_c, direction)
			second_branch = try_short_tpl_ext(hash_table, second_query, direction, 5)

			// If there is probable branching
			if (main_branch != second_branch): 
				second_branch = try_short_tpl_ext(hash_table, second_query, direction, hash_table.read_length - 5)
				has_junction_read = find_junction_reads(template, second_branch, hash_table, direction)
				// If the second branch is supported by some junction reads, extend it
				if (has_junction_read):
					print "Extending the main template first"
					branching_locus = template.length
					// Not using the result of try_short_tpl_ext, because their kmers are not marked as used.
					shift_template(template, c, direction)
					// Recursivly call this function itself to extend
					template_ext(template, hash_table, main_query, direction)
					print "Extending the branch at %d" % branching_locus
					// For the second branch, extend to only one direction, because another direction is extended by main template already
					// Not using the result of try_short_tpl_ext, because their kmers are not marked as used.
					second_branch = second_c
					// In case of branching on short templates
					second_branch.virtual_tail = cut_template_tail(template, hash_table.read_length - 5, direction)
					template_ext(second_branch, hash_table, second_query, direction)
					branching_points.add((template, second_branch, branching_locus, direction))
					ALL_TEMPLATES.add(second_branch)
					// Jump out of the function
					return
			else:
				// It is just an ambiguous base because of sequencing errors, etc. Record it, continue extension
				mark_kmer_used(hash_table, branch_query)
				AMBIGUOUS_BASES.add((template, template.length))

		// For every kmer used during this extension, mark the frequency as 0.
		mark_kmer_used(hash_table, query)
		// Extend template and the query by one char
		shift_template(template, c, direction)
		query = shift_query(query, c, direction)

/**
	Quickly extend a query to some length (parameter max_len), does not mark kmers as used
	Return: the short branch template
**/
Function try_short_tpl_ext(hash_table, query, direction, max_len):
	branch_template = 0
	query_copy = query
	for (i = 0; i < max_len; i++):
		counters = count_next_kmers(hash_table, query_copy, direction)
		c = get_max_index(counters)
		if (c == -1):
			break;
		else:
			shift_template(branch_template, c, direction)
			query_copy = shift_query(query_copy, c, direction)
	return branch_template	

/**
	Check whether there are reads in the junction area
	Return: 1 if there is at least one junction read, 0 otherwise
	E.g.:	read_length: 15
			short_shift: 5 (hard coded now)
			main_template:     ATC...CCCCCCCCCCCCCCC (length 1000)
			branch:                                 TTTTTTTTTT (length 15-5=10)
			junction_template:            CCCCCCCCCCTTTTTTTTTT (length 15*2-5*2=20)
			junction_reads:                CCCCCCCCCTTTTTT
			                                  CCCCCCTTTTTTTTT 
**/
Function find_junction_reads(main_template, branch, hash_table, direction):
	// Get the tail of the main template
	main_template_tail = cut_template_tail(main_template, hash_table.read_length - 5, direction)
	// Concat the two branches, temporarily
	probable_junction_template = concat(main_template_tail, branch)
	junction_reads = []
	for (i = 0; i < probable_junction_template.length - hash_table.read_length; i++):
		// query length is same as read length, allowing 1 mismatch
		query = get_query_int(probable_junction_template, i, hash_table.read_length)
		junction_reads.add(align(hash_table, query, 1))
	if (junction_reads.length > 0):
		return 1
	return 0


// ========= Supportive functions ==========

/**
	Get the tail of a template. 
	If the length is not long enough, try to get from the 'virtual_tail'
**/
Function cut_template_tail(template, tail_length, direction):
	tail = NULL
	tpl_copy = template
	if (tpl_copy.length < tail_length):
		tpl_copy = concat(template.virtual_tail, tpl_copy)
	return simple_tail(tpl_copy, tail_length, direction)

/**
	Start from locus i on the template, return the int64 of the new query with certain length
	E.g.:	template: ATC...TTTTTCCCCC (length 1000)
			i: 992
			length: 5
			to return:        TTTCC (11 11 11 01 01)
**/
Function get_query_int(template, i, length):
	// Return the int64 value of the new query

/** 
	Simply get index of second largest value
**/
Function get_second_freq_char(counters, max):
	counters_without_max = []
	for (i = 0; i < 4; i++):
		if (i == max):
			counters_without_max[i] = 0
		else:
			counters_without_max[i] = counters[i]
	return get_max_index(counters_without_max)

/**
	In an array, if some value is maximum, return its index
**/
Function get_max_index(counters):
	for (i = 0; i < 4; i++):
		if (counters[i] > max):
			max = counters[i]
			max_index = i
	// If no probable next char, return -1
	if (max == 0)
		return -1
	return max_index


// ========= Operations on the hash table ==========

/**
	Get the frequencies of all probable next kmers
	Return: [frequencies of next kmers]
	E.g.:	query: 			AATCT
			frequencies: 	 ATCTA: 10
							 ATCTC: 2
							 ATCTG: 0
							 ATCTT: 1
			to return: 		[10, 2, 0, 1]
**/
Function count_next_kmers(hash_table, query, direction):
	counters = []
	next_probable_kmer = 0									
	for (i = 0; i < 4; i++):
		counters[i] = 0
		next_probable_kmer = shift_bit(query, i, direction)
		counters[i] = get_kmer_count(hash_table, next_probable_kmer)
		// Get its reverse complement as well
		next_probable_kmer = rev_com_kmer(next_probable_kmer)
		counters[i] += get_kmer_count(hash_table, next_probable_kmer)
	return counters

Function align(hash_table, query, n_mismatch):
	// Return a list of reads containing the query, allowing n_mismatch

Function mark_kmer_used(hash_table, query):
	hash_table[query] = NULL
	hash_table[rev_com_kmer(query)] = NULL
