#itreecount

itreecount counts reads from bam files on genomic regions supplied trough a GTF
file. The processing is parallelized over chromosomes in the bam file only. Counting 
is done identical to the union method of the popular (but slow) htseq-count python package.

##Dependencies
itreecount depends on the following perl packages:

	Parallel::ForkManager  (CPAN)
	Bio::DB::Sam           (CPAN)
	Set::IntervalTree      (CPAN)
	Getopt::Long           (core)
	File::Temp             (core)

##Usage
    itreecountpe.pl --gtf <gtffile> --ncpu 10 --output counts.txt <bamfile.bam>

The bam file can be single or paired-end, needs to be indexed (and not namesorted).

##Ouput

Tab delimited text file with gene identifiers from the GTF file followed by the read counts. Also per chromosome the ambigious, not unique and no feature alignments are reported.

##Warnings:

The read pairs in the bam file are counted once (if they fall on the same gene), but because the program parses position sorted data it temporarily stores these reads in memory. If the dataset has many mated reads that are separated over a large distance the memory consumption might grow.

