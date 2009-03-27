#!/usr/bin/perl -w

use FileHandle;

my $minid = 0;
my $maxid = 100;
my $dupid = 98;

my $usage = "Usage: $0 <Stockholm database> <output file>\n";
$usage .= "            [-minid <minimum pairwise \%ID>]        (default is $minid)\n";
$usage .= "            [-maxid <maximum pairwise \%ID>]        (default is $maxid)\n";
$usage .= "            [-dupid <maximum nonredundant \%ID>]    (default is $dupid)\n";
$usage .= "\n";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-h") {
	die $usage;
    } elsif ($arg eq "-minid") {
	defined ($minid = shift) or die $usage;
    } elsif ($arg eq "-maxid") {
	defined ($maxid = shift) or die $usage;
    } elsif ($arg eq "-dupid") {
	defined ($dupid = shift) or die $usage;
    } else {
	push @argv, $arg;
    }
}
die $usage unless @argv == 2;
my ($stockfile, $outfile) = @argv;

my (@stock, $stock);

warn "Reading alignments from $stockfile\n";
my $in = FileHandle->new;
$in->open ("<$stockfile") or die "Couldn't open $stockfile: $!";

warn "Writing alignments to $outfile\n";
my $out = FileHandle->new;
$out->open (">$outfile") or die "Couldn't open $outfile: $!";
$out->autoflush (1);

local $_;
while ($_ = $in->getline()) {
    next unless /\S/;
    if (/^\s*(\S+)\s+(\S+)\s*$/) {  # sequence
	my ($name, $seq) = ($1, $2);
	$stock = addTag (\@stock, $stock, $name, "", $seq);
    } elsif (/^\s*\#=GR\s*(\S+)\s*(\S+)\s*(\S+)/) {  # by-seq by-col (#=GR)
	my ($name, $tag, $data) = ($1, $2, $3);
	$stock = addTag (\@stock, $stock, $name, $tag, $data);
    } elsif (/^\s*\#=GC\s*(\S+)\s*(\S+)/) {  # by-col (#=GC)
	my ($tag, $data) = ($1, $2);
	$stock = addTag (\@stock, $stock, "", $tag, $data);
    } elsif (/\s*\/\/\s*$/) {
	splitStock ($stock, $out);
	$stock = undef;
    }
}
splitStock ($stock, $out) if defined $stock;

$out->close() or die "Couldn't close $outfile: $!";
$in->close();

sub getPair {
    my ($stock, $i, $j, $seqs, $quick) = @_;
    my ($iname, $jname) = @$seqs[$i,$j];
    my ($itagval, $jtagval) = ($stock->{$iname}, $stock->{$jname});
    my ($igapped, $jgapped) = (\$itagval->{""}, \$jtagval->{""});
    if (length($$igapped) ne length($$jgapped)) {
	warn "Row data for ($iname,$jname) have different lengths (", length($$igapped), ",", length($$jgapped), "); skipping\n";
	return (undef, undef, 0, {});
    }
    my ($iseq, $jseq) = ("", "");
    my @itags = grep (length(), keys %$itagval);
    my @jtags = grep (length(), keys %$jtagval);
    my %itag = map (($_ => ""), @itags);
    my %jtag = map (($_ => ""), @jtags);
    my $gctags = exists ($stock->{""}) ? $stock->{""} : {};
    my @gctags = keys %$gctags;
    my %pairgctags = map (($_ => ""), @gctags);
    my ($match, $total) = (0,0);
    for (my $col = 0; $col < length($$igapped); ++$col) {
	my $ic = uc substr($$igapped,$col,1);
	my $jc = uc substr($$jgapped,$col,1);
	my $igap = isGap ($ic);
	my $jgap = isGap ($jc);
	my $colEmpty = $igap && $jgap;
	my %gccol;
	foreach my $gctag (@gctags) {
	    my $gcval = $$gctags{$gctag};
	    my $gccol = substr ($gcval, $col, 1);
	    $gccol{$gctag} = $gccol;
	    if (!isGap ($gccol)) {
		$colEmpty = 0;
	    }
	}
	unless ($colEmpty) {
	    if (!$igap && !$jgap) {
		if ($ic eq $jc) { ++$match }
		++$total;
	    }
	    unless ($quick) {
		$iseq .= $ic;
		$jseq .= $jc;
		foreach my $itag (@itags) { $itag{$itag} .= substr ($itagval->{$itag}, $col, 1) }
		foreach my $jtag (@jtags) { $jtag{$jtag} .= substr ($jtagval->{$jtag}, $col, 1) }
		foreach my $gctag (@gctags) { $pairgctags{$gctag} .= $gccol{$gctag} }
	    }
	}
    }
    my $id = $total ? 100 * $match / $total : 0;
    return ([$iname,$iseq,\%itag], [$jname,$jseq,\%jtag], $id, \%pairgctags);
}

sub splitStock {
    my ($stock, $out) = @_;
    my @seqs = grep (length, keys %$stock);
    my $nseqs = @seqs;
    return if $nseqs == 0;
    warn " [splitting a ", @seqs+0, "-row alignment]\n";
    my @index;
    my @max;
    for (my $i = 0; $i < @seqs; ++$i) {
#	my $i = $seqs[$a];
	my $discard = 0;
	for (my $b = $i+1; $b < @seqs; $b++) {
	    #my $j = $b;
	    my ($iref, $jref, $id, $gctags) = getPair ($stock, $i, $b, \@seqs, 1); 
	    if ($id >= $dupid) { $discard = 1; last } #removed || $id < $minid test 4/10/2008 LEB
	    $max[$i] = $id if ((!defined($max[$i])) || ($id > $max[$i]));  #new max test, 4/14/2008 LEB
	    $max[$b] = $id if ((!defined($max[$b])) || ($id > $max[$b]));  # ..
	}
	if ($max[$i] < $minid) { $discard = 1;}
	push @index, $i unless $discard;
	warn " [processed ", $i+1, "/$nseqs sequences, kept ", @index+0, "; last max \%ID $max[$i]\]\n" unless $discard;
    }
    warn " [trimmed alignment from ", @seqs+0, " to ", @index+0, " rows]\n";

    if (!@index) {warn " [no alignments to process, skipping]\n"; return}
    my $weight = 1 / @index;
    for (my $a = 0; $a < @index; ++$a) {
	my $i = $index[$a];
	my @output;
	for (my $b = $a+1; $b < @index; ++$b) {
	    my $j = $index[$b];
	    my ($iref, $jref, $id, $gctags) = getPair ($stock, $i, $j, \@seqs, 0);
	    my ($iname, $iseq, $itags) = @$iref;
	    my ($jname, $jseq, $jtags) = @$jref;
	    if ($id >= $minid && $id <= $maxid) {
		
		printf STDERR " [pairwise alignment ($iname,$jname) has weight %.2g, identity %d%%]\n", $weight, int($id+.5);
		
		while (length($iname)<length($jname)) { $iname .= " " }
		while (length($jname)<length($iname)) { $jname .= " " }

		push @output, "#=GF WT $weight\n";
		push @output, "$iname $iseq\n";
		push @output, "$jname $jseq\n";
		foreach my $itag (keys %$itags) { push @output, "#=GR $iname $itag $$itags{$itag}\n" }
		foreach my $jtag (keys %$jtags) { push @output, "#=GR $jname $jtag $$jtags{$jtag}\n" }
		foreach my $gctag (keys %$gctags) {
		    my $gcval = $$gctags{$gctag};
		    my $gctag_hdr = "#=GC $gctag";
		    while (length ($gctag_hdr) < length($iname)) { $gctag_hdr .= " " }
		    push @output, "$gctag_hdr $gcval\n";
		}
		push @output, "//\n";
	    }
	}
	$out->print (@output);
    }
}

sub addTag {
    my ($stockArray, $stock, $name, $tag, $data) = @_;
    if (!defined ($stock)) {
	$stock = {};
	push @$stockArray, $stock;
    }
    if (!exists $stock->{$name}) {
	$stock->{$name} = {};
    }
    if (!exists $stock->{$name}->{$tag}) {
	$stock->{$name}->{$tag} = "";
    }
    $stock->{$name}->{$tag} .= $data;
    return $stock;
}

sub isGap { my $c = shift; return $c eq '.' || $c eq '-' }
