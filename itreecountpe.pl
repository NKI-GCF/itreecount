#!/usr/bin/perl 
#===============================================================================
#
#         FILE: itreecount.pl
#
#        USAGE: ./itreecount.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Arno Velds (), a.velds@nki.nl
#      COMPANY: NKI
#      VERSION: 1.0
#      CREATED: 02/15/2013 10:26:14 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;

use Bio::DB::Sam;
use Set::IntervalTree;
use Getopt::Long;
use Data::Dumper;

my $gtffile;
my $chr;

GetOptions("gtf:s"=>\$gtffile, "chr:s"=>\$chr) or die "Usage";
my $bamfile = shift;

(my $chrtree, my $genelist) = read_gtf($gtffile, $chr);
#prepare the count structure
my %counts = map {$_ => 0} keys $genelist;
#the wait structure
my %delayed = ();

print STDERR scalar(keys(%$genelist)), " genes found\n";

my $bam          = Bio::DB::Bam->open($bamfile);
my $header       = $bam->header;
my $target_count = $header->n_targets;
my $target_names = $header->target_name;

my $index = Bio::DB::Bam->index_open_in_safewd($bamfile) or die "Indexed bams only\n";

my $callback = sub {
	my $alignment = shift;
	my $start       = $alignment->start;
	my $end         = $alignment->end;
	my $seqid       = $target_names->[$alignment->tid];
	#check for pairs?
	
	my $flag = $alignment->flag;
	#we don't get unmapped reads in tophat output, but skip em anyway
	next if $flag & 4;

	my @countthis;

	if ($flag & 1) {
		#paired end
		if ($flag & 8) {
			#mate is unmapped treat this read as single end
			push @countthis, $alignment;

		} else { #mate is mapped
			#is the mate on the same chromosome? if not than this read is ambiguous
			++$counts{ambiguous} && return if $alignment->mtid != $alignment->tid;
			
			#do we already have the mate info?
			my $mateid = join(":", $alignment->qname, $alignment->mate_start);
			#print STDERR "$mateid\n";
			if (exists $delayed{$mateid}) {
				my $a2 = shift(@{$delayed{$mateid}});
				push @countthis, ($alignment, $a2);
				delete $delayed{$mateid} if $#{$delayed{$mateid}} == -1;
			} else {
				#delay this read until we see the mate
				push @{$delayed{join(":", $alignment->qname, $start)}}, $alignment;
				return;
			}
		}
	} else {
		#unpaired read. Just count it normally
		push @countthis, $alignment;
	}


	my %genes;
	foreach my $a (@countthis) {
		#multimapper (one or both...don't really care
		++$counts{alignment_not_unique} && return if $a->aux_get("NH") > 1;

		#get the cigar string
		my $cigarray = $a->cigar_array;
		my $pos = $a->start;
		foreach my $e (@$cigarray) {
			if ($e->[0] =~ "^M") {
				$genes{$_}++ foreach (@{$chrtree->fetch($pos, $pos + $e->[1]);});
				$pos += $e->[1];
			} elsif ($e->[0] =~ "^N") {
				$pos += $e->[1];
			} elsif ($e->[0] =~ "^I") {
				#insertion in reference (don't count)
			} elsif ($e->[0] =~ "^D") {
				#deletion in reference (count)
				$pos += $e->[1];
			} else {
				warn "Cigar operation $e->[0] not supported", join("", map {$_->[0] . $_->[1]} @$cigarray), "\n";
			}
		}
	}

	my @ugenes = keys %genes;
	++$counts{no_feature} && return if $#ugenes == -1;
	++$counts{ambiguous} && return if $#ugenes > 0;
	$counts{$ugenes[0]}++;
	return;
};

$index->fetch($bam,$header->parse_region($chr),$callback);

print STDERR scalar(keys(%delayed)), " mates never matched\n";
$counts{alignment_not_unique} += scalar(keys(%delayed));

print join("\t", $_, $counts{$_}),"\n" foreach sort keys %counts;
exit;


sub read_gtf {
	my $file = shift;
	my $chr = "" . shift;

	open(IN,"<", $file) or die "$!\n";

	my $tree = Set::IntervalTree->new;
	my $haveread = 0;
	my %genes;
	my $nread = 0;
	while(my $line = <IN>) {
		my @e = split /\t/, $line;
		if ($e[0] ne $chr) {
			last if($haveread);
			next;
		}
		next unless $e[2] eq "exon";
		next unless $e[8] =~ m/gene_id\s\"(\w+)\"/;
		my $geneid = $1;
		$genes{$geneid} = 0;
		$tree->insert($geneid, $e[3], $e[4]+1);
		
		$haveread = 1;
		$nread++;
	}
	print STDERR "Read $nread lines from GTF\n";
	return ($tree, \%genes);
}




