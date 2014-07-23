#!/usr/bin/env perl 

use strict;
use warnings;

use Parallel::ForkManager;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
use Cwd 'abs_path';

use FindBin;
use lib "$FindBin::Bin";
use TreeCount;

my $gtffile;
my $ncpu = 10;
my $verbose;
my $output;

GetOptions("gtf:s"=>\$gtffile, "ncpu|t:i"=>\$ncpu, "verbose|v"=>\$verbose, "output|o:s"=>\$output) or die "Usage";
my $bamfile = abs_path(shift);

(my $iforest, my $genelist) = read_gtf_as_intervalforest($gtffile);

my ($bam, $header) =  open_bam($bamfile);
my $target_count = $header->n_targets;
my $target_names = $header->target_name;
my $chrlengths = $header->target_len;

#prepare the fork manager
my $pm = Parallel::ForkManager->new($ncpu);

#the location for the temp files
my $tmpdir = tempdir( CLEANUP => 1 );
my %tmpfiles;

#process the chromosomes in large to small order should give best runtime
my @o = sort { $chrlengths->[$b] <=> $chrlengths->[$a] } 0 .. ($target_count - 1);
foreach my $tid (@o) {
	my $chr =  $target_names->[$tid];
	print STDERR "Forking child for chromosome $chr\n";

	#open the temp file for this chromosome
	my ($fh, $filename) = tempfile( UNLINK => 0, DIR => $tmpdir );

	if($verbose) {
		print STDERR "Child $$ writing temporary data to $filename\n";
	}

	$tmpfiles{$target_names->[$tid]} = $filename;

	$pm->start and next; # do the fork

	#we must reopen the bam in the forked child
	my ($fbam, $fheader, $findex) = open_bam($bamfile);

	#initialize all counts at zero
	my $count = { map {$_ => 0} keys %{$genelist->{$chr}} };
	my $delayed = {};

	#load the requested chromsome from the bam file and call the callback for every read
	$findex->fetch($fbam,$fheader->parse_region($chr), \&TreeCount::count_read_callback, [$iforest, $count, $delayed, $target_names]);

	print STDERR "Chromosome $chr done. ", scalar(keys(%$delayed)), " mates never matched\n" if $verbose;
	print $fh join("\t", $_, $count->{$_}),"\n" foreach sort keys %$count;
	close($fh);

	$pm->finish; # do the exit in the child process
}
$pm->wait_all_children;

#write report
print STDERR "Collecting data parts\n";
if(defined $output) {
	open(OUT,">", $output) or die "$!\n";
	select OUT;
}

#collect all data
for my $ch (@$target_names) {
	#get the formatted path
	open(I,"<", $tmpfiles{$ch}) or warn "Temp file $tmpfiles{$ch} for chromosome $ch not found?\n";
	print while <I>;
	close(I);
}

exit;


