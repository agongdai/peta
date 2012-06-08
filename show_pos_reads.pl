#!/usr/bin/perl
#use warnings;
use 5.010;
use File::Basename;

$rnaseq_file = "read/rnaseq_200.fa";
if (@ARGV == 3) {
	$rnaseq_file = $ARGV[0];
	$tx_id = $ARGV[1];
	$pos = $ARGV[2];
} else {
	if (@ARGV == 2) {
		$tx_id = $ARGV[0];
		$pos = $ARGV[1];
	} else {
		say "Format: ./show_gene_no.pl [rnaseq file name] read_id pos";
		exit;
	}
}

open RNASEQ,    	"<", $rnaseq_file;
say $pos;

$counter = 0;
$print_seq = 0;
$range_low = $pos - 60;
$range_high = $pos + 60;
say "======================================";
say "Reads nearby: ";

open RNASEQ,    	"<", $rnaseq_file;
while (<RNASEQ>) {
	chomp;
	if (/>/) {
		@fields = split("\t");
		if ($fields[1] eq $tx_id && $fields[2] > $range_low && $fields[2] < $range_high) {
			$print_seq = 1;
			print;
			if ($fields[2] < $pos) {
				print("========");
			}
			print("\n");
		}
	} else {
		if ($print_seq) {
			say;
			$print_seq = 0;
		}
	}
}
close(RNASEQ);
