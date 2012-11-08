#!/usr/bin/perl
#use warnings;
use 5.010;
use File::Basename;

open GENES_LIST,    	"<", "read/gene_10.txt";
open GENES_EDGES, 		">", "graph/gene_10_edges.txt";
open GENES_LEN, 		">", "graph/gene_10_lengths.txt";

%genes;				# Hash map: "<start_pos>,<end_pos> => <node_id>"
$start_pos = 123;
$end_pos = 123;
$gene = "";
$node_id = 0;
$i = 0;
$gene_nodes = "";
$totle_len = 0;
$len = 0;
@starts;
@ends;
@stand_alone;

while (<GENES_LIST>) {
	chomp;
	if (/#/) {
	} else {
		$gene_nodes = "";
		$totle_len = $i = 0;
		@fields = split("\t");
		# say "Dealing with gene " . $fields[1];
		$start_str = $fields[9];
		$end_str = $fields[10];
		@starts = split(",", $start_str);
		@ends = split(",", $end_str);
		foreach(@starts) {
			$gene = $_ . "," . $ends[$i];
			$len = $ends[$i] - $_;
			$totle_len += $len;
			if (!$genes{$gene}) {
				$genes{$gene} = $node_id++;
			}
			say GENES_LEN $genes{$gene} . ":$len:[$_," . $ends[$i] . "]";
			$gene_nodes .= $genes{$gene} . "[$len],";
			$i++;
		}
		chop($gene_nodes);
		say GENES_EDGES "#" . $fields[1] . "|len=$totle_len";
		say GENES_EDGES $gene_nodes;
	}
}

close(GENES_LIST);
close(GENES_EDGES);
close(GENES_LEN);
