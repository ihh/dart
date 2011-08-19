package SequenceIterator;

# Ian Holmes  --  ihh@sanger.ac.uk  --  Aug 26 1998
# Updated 2009, 2011 -- ihh@berkeley.edu

#
# Perl OOP includes
#

use vars qw($AUTOLOAD @ISA @EXPORT_OK %EXPORT_TAGS);
use base qw(Exporter);
%EXPORT_TAGS = ('all' => [@EXPORT_OK = qw(itertemp iterseq biterseq gbiterseq printseq revcomp)]);

#
# package variables
#

($SequenceIterator::prog = $0) =~ s#.*/([^/]+)$#$1#;
$SequenceIterator::tempPrefix = "$ENV{HOME}/tmp/$SequenceIterator::prog.$$.";
$SequenceIterator::tempCount = 0;

#
# cleanup routine
#

END {
    foreach (keys %SequenceIterator::tempFile) {
	if (-e $_) { unlink $_ or warn "Couldn't unlink $_: $!" }
    }
}
$SIG{INT} = $SIG{KILL} = sub { die "\n" };

#
# itertemp - subroutine to iterate through every sequence in a FASTA database, making temporary files
#

sub itertemp {
    my ($file,$sub,$filter) = @_;

    local *file;
    local *tempFile;

    my $tempFile = $SequenceIterator::tempPrefix . ++$SequenceIterator::tempCount;
    $SequenceIterator::tempFile{$tempFile}++;

    local *FILE;
    if (defined $filter) {
	open FILE, "cat $file | $filter |" or die "$0: couldn't open $file: $!";
    } else {
	open FILE, $file or die "$0: couldn't open $file: $!";
    }

    my $sep = $/;
    $/ = ">";
    my $dummy = <FILE>;
    my $name;
    while (($/ = "\n", $name = <FILE>)[1]) {
	chop $name;
	$/ = ">";
	my $seq = <FILE>;
	chop $seq;
	open tempFile, ">$tempFile" or die "$0: couldn't open $tempFile for writing: $!";
	print tempFile ">$name\n$seq";
	close tempFile;
	$/ = $sep;
	$seq =~ s/\s//g;
	&$sub($tempFile,$name,$seq);
	unlink $tempFile;
    }
    close FILE;
    $/ = $sep;
    delete $SequenceIterator::tempFile{$tempFile};
}


#
# iterseq - same as itertemp, but doesn't create a temporary file
#

sub iterseq {
    my ($file,$sub,$filter) = @_;

    local *FILE;
    if (defined $filter) {
	open FILE, "cat $file | $filter |" or die "$0: couldn't open $file: $!";
    } else {
	open FILE, "<$file" or die "$0: couldn't open $file: $!";
    }

    my $sep = $/;
    $/ = ">";
    my $dummy = <FILE>;
    my $name;
    while (($/ = "\n", $name = <FILE>)[1]) {
	chop $name;
	$/ = ">";
	my $seq = <FILE>;
	$seq = "" unless defined $seq;
	chomp $seq;
	$/ = $sep;
	$seq =~ s/\s//g;
	&$sub($name,$seq);
    }
    close FILE;
    $/ = $sep;
}

#
# biterseq - Bioperl iterseq
#

sub biterseq {
    my ($file, $format, $sub) = @_;
    my $in = Bio::SeqIO->new ('-file' => $file, '-format' => $format);
    while (my $seq = $in->next_seq) { &$sub ($seq) }
}

#
# gbiterseq - GenBank biterseq
#

sub gbiterseq { biterseq (shift, "GenBank", shift) }

#
# printseq - print a sequence
#

sub printseq {
    my $name = shift;
    my $seq = shift;
    my $fh = @_ ? shift : \*STDOUT;
    my $width = @_ ? shift : 50;
    my $i;
    print $fh ">$name\n";
    for ($i=0;$i<length $seq;$i+=$width) { print $fh substr($seq,$i,$width),"\n" }
}

# revcomp
sub revcomp {
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/acgtuACGTU/tgcaaTGCAA/;
    return $seq;
}

1;


__END__

=head1 NAME

SequenceIterator.pm - External iterator for FASTA databases

=head1 DESCRIPTION

B<SequenceIterator.pm>
This module provides a subroutine called 'itertemp',
that takes as arguments a FASTA database filename
and a subroutine reference.
For every sequence in the database,
a temporary file is created containing that sequence alone
and the supplied subroutine is called with the
filename of the temporary file as its first argument.
This is a classical external iterator.

An alternative formation is 'iterseq'.
This is the same as 'itertemp', but creates no temporary file.

=head2 Example

A skeletal implementation of 'blastdb':

    #!/usr/bin/perl
    
    use SequenceIterator qw(all);

    ($executable,$database,$querydb) = @ARGV;

    itertemp ($querydb, sub {
                              my $query = shift;
                              system "$executable $database $query";
  			   });
        

=head1 METHOD CALLS

=over 5

=item iterseq

iterseq($database,$sub)

iterseq($database,$sub,$filter)

Calls &$sub($sequence_name,$sequence_data) for
every sequence in $database.

If a filter program is specified in $filter then $database
is piped through that first.

=item itertemp

itertemp($database,$sub)

itertemp($database,$sub,$filter)

Creates a temporary file for each successive sequence in $database
and then calls &$sub($temp_file_name,$sequence_name,$sequence_data).

If a filter program is specified in $filter then $database
is piped through that first.

=item revcomp

my $rev = revcomp($seq)

Reverse-complements a sequence.


