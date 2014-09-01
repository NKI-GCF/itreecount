package TreeCount;

use strict;
use warnings;

use Bio::DB::Sam;
use Set::IntervalTree;
use Data::Dumper;

our (@ISA, @EXPORT);
BEGIN {
	require Exporter;
	@ISA = qw(Exporter);
	@EXPORT = qw(open_bam read_gtf_as_intervalforest);
}

#name the constants
my $AMBIGUOUS = "ambiguous";
my $ALIGNMENT_NOT_UNIQUE = "alignment_not_unique";
my $NO_FEATURE = "no_feature";

sub open_bam  {
	my $file = shift;

	my $bam          = Bio::DB::Bam->open($file);
	my $header       = $bam->header;
	my $target_count = $header->n_targets;
	my $target_names = $header->target_name;
	
	my $index = Bio::DB::Bam->index_open_in_safewd($file) or die "Indexed bams only, failed on $file\n";

	return ($bam, $header, $index);
}

sub count_read_callback {
	my $alignment = shift;

	my $data = shift;
	my ($iforest, $counts, $delayed, $chr_names, $stranded) = @$data;

	my $start = $alignment->start;
	my $end = $alignment->end;
	my $chr = $chr_names->[$alignment->tid];
	
	my $flag = $alignment->flag;
	#we don't get unmapped reads in tophat output, but skip them anyway
	next if $flag & 4;

	#also skip non primary alignments, vendor-failed and marked (optical) duplicates
	next if ($flag & 256) && ($flag & 512) && ($flag & 1024);

	my @countthis;

	if ($flag & 1) {
		#paired end
		if ($flag & 8) {
			#mate is unmapped treat this read as single end
			push @countthis, $alignment;
		} else { #mate is mapped
			#is the mate on the same chromosome? if not than this read is ambiguous
			++$counts->{"$AMBIGUOUS-$chr"} && return if $alignment->mtid != $alignment->tid;
			
			#do we already have the mate info?
			my $mateid = join(":", $alignment->qname, $alignment->mate_start);
			if (exists $delayed->{$mateid}) {
				my $a2 = shift(@{$delayed->{$mateid}});
				push @countthis, ($alignment, $a2);
				delete $delayed->{$mateid} if $#{$delayed->{$mateid}} == -1;
			} else {
				#delay this read until we see the mate
				push @{$delayed->{join(":", $alignment->qname, $start)}}, $alignment;
				return;
			}
		}
	} else {
		#unpaired read. Just count it normally
		push @countthis, $alignment;
	}


	my %genes;
	foreach my $a (@countthis) {
		#multimapper (one or both...don't really care)  (for non tophat results cut on mapq)
		if (defined $a->aux_get("NH") ) {
			++$counts->{"$ALIGNMENT_NOT_UNIQUE-$chr"} && return if $a->aux_get("NH") > 1
		} else {
			++$counts->{"$ALIGNMENT_NOT_UNIQUE-$chr"} && return if $a->qual < 10
		}

		#get the cigar string
		my $cigarray = $a->cigar_array;
		#the cigar parser parses only a subset of the possible elements. This
		#should be enough to process TopHat data, but might need enhancements when
		#other aligners are used.
		my $pos = $a->start;
		foreach my $e (@$cigarray) {
			if ($e->[0] =~ /^M/) {
				#store gene name if segment overlaps a gene (optionally stranded)
				foreach ((@{$iforest->{$chr}->fetch($pos, $pos + $e->[1]);})) {
					if($_->[0] eq "ENSG00000008735") {
					$genes{$_->[0]}++ if !$stranded || ( (($flag & 16) != 0 ) == $_->[1]);
				}
				$pos += $e->[1];
			} elsif ($e->[0] =~ /^N/) {
				$pos += $e->[1];
			} elsif ($e->[0] =~ /^I/) {
				#insertion in reference (don't count)
			} elsif ($e->[0] =~ /^D/) {
				#deletion in reference (count)
				$pos += $e->[1];
			} elsif ($e->[0] =~ /^S/ || $e->[0] =~ /^H/) {
				#hard/soft clipped based we ignore
			} else {
				warn "Cigar operation $e->[0] not supported", join("", map {$_->[0] . $_->[1]} @$cigarray), "\n";
			}
		}
	}

	my @ugenes = keys %genes;
	++$counts->{"$NO_FEATURE-$chr"} && return if $#ugenes == -1;
	++$counts->{"$AMBIGUOUS-$chr"} && return if $#ugenes > 0;
	$counts->{$ugenes[0]}++;
	return;
}

sub read_gtf_as_intervalforest {
	my $file = shift;

	print STDERR "Parsing GTF file\n";
	open(IN,"<", $file) or die "$!\n";

	my %forest;
	my %genes;
	my $nread = 0;
	my %genecount;
	while(my $line = <IN>) {
		my @e = split /\t/, $line;
		
		next unless $e[2] eq "exon";
		next unless $e[8] =~ m/gene_id\s\"(\w+)\"/;

		my $chr = $e[0];
		if (!exists $forest{$chr}) {
			$forest{$chr} = Set::IntervalTree->new;
		}
		my $geneid = $1;
		$genes{$chr}->{$geneid} = 0;
		$genecount{$geneid} = 0;
		$forest{$chr}->insert([$geneid, $e[6] eq "+"], $e[3], $e[4]+1);
		
		$nread++;
	}
	close(IN);
	print STDERR "Read $nread exons in ", scalar(keys(%genecount)), " genes from GTF\n";
	return (\%forest, \%genes);
}



1;

