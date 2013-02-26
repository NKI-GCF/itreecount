use strict;
use warnings;


sub read_gtf_as_intervalforest {
	my $file = shift;
	my $chr = "" . shift;

	open(IN,"<", $file) or die "$!\n";

	my %forest;
	my $tree = Set::IntervalTree->new;
	my %genes;
	my $nread = 0;
	while(my $line = <IN>) {
		my @e = split /\t/, $line;
		if (!exists $forest{$chr}) {
			$forest{$chr} = Set::IntervalTree->new;
		}
		next unless $e[2] eq "exon";
		next unless $e[8] =~ m/gene_id\s\"(\w+)\"/;
		my $geneid = $1;
		$genes{$geneid} = 0;
		$forest{$chr}->insert($geneid, $e[3], $e[4]+1);
		
		$nread++;
	}
	print STDERR "Read $nread lines from GTF\n";
	return (\%forest, \%genes);
}



1;

