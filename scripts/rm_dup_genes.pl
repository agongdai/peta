#!/usr/bin/perl
#use warnings;
use 5.010;
use File::Basename;

open ORI_GENES,    	"<", "read/gene_10_ori.txt";
open GENES,    		">", "read/gene_10.txt";

%genes;
say "Duplicated genes removed: ";

while (<ORI_GENES>) {
	chomp;
	if (/#/) {
		say GENES;
	} else {
		@fields = split("\t");
		$gene_name = $fields[1];
		if (!$genes{$gene_name}) {
			if (/MIR/) {
				say;
			} else {
				$genes{$gene_name} = $_;
				say GENES;
			}
		} else {
			say;
		}
	}
}
close(ORI_GENES);
close(GENES);
