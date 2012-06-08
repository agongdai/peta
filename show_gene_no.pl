#!/usr/bin/perl
#use warnings;
use 5.010;
use File::Basename;

open READS, 			">", "read/start_reads_gene.txt";
open SEQ, 			">", "read/start_reads_seq.txt";

$rnaseq_file = "read/rnaseq_200.fa";
if (@ARGV >= 2) {
	# $rnaseq_file = $ARGV[0];
} else {
	if (@ARGV == 1) {
		$read_id = $ARGV[0];
	} else {
		say "Format: ./show_gene_no.pl [rnaseq file name] read_id";
		exit;
	}
}

$counter = 0;
%target_ids;
foreach $id (@ARGV) {
	$target_ids{$id} = 1;
}

open RNASEQ,    	"<", $rnaseq_file;
while (<RNASEQ>) {
	chomp;
	if (/>/) {
		$print_seq = 0;
		if ($target_ids{$counter} == 1) {
			@fields = split("\t");
			$id = ">" . $counter;
			if ($fields[0] eq $id) {
				say;
				say SEQ $_;
				$tx_id = $fields[1];
				say READS $tx_id;
				$print_seq = 1;
			}
		}
		$counter++;
	} else {
		if ($print_seq) {
			say;
			say SEQ $_;
		}
	}
}
close(RNASEQ);
close(READS);
close(SEQ);
