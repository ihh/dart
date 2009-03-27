#!/usr/bin/perl -w

use Stockholm;

my $amino = 0;
my $slim = 0;

my $usage = "Usage: $0 [-amino] [-slim] [<files...>]";

my @argv;
while (@ARGV) {
    my $arg = shift;
    if ($arg eq "-amino") {
	$amino = 1;
    } elsif ($arg eq "-slim") {
	$slim = 1;
    } elsif ($arg =~ /^-.+/) {
	die $usage;
    } else {
	push @ARGV, $arg;
    }
}

my @pandit;
while (<>) {
    chomp;
    if (/^\s*\/\/\s*$/) {
	print_stock (\@pandit);
	@pandit = ();
    } else {
	push @pandit, $_;
    }
}
print_stock (\@pandit) if @pandit;

# convert an array of PANDIT lines into Stockholm
sub print_stock {
    my ($pandit) = @_;

    my $stock = Stockholm->new;
    my ($name, $aln, $dln);
    foreach my $line (@$pandit) {
	if ($line =~ /^(\S\S\S)\s+(.*)$/) {
	    my ($field, $data) = ($1, $2);
	    $field = uc $field;

	    if ($field eq 'FAM') {
		push @{$stock->gf_AC}, $data;
		$aln = $dln = undef;

	    } elsif ($field eq 'PID') {
		push @{$stock->gf_ID}, $data;

	    } elsif ($field eq 'DES') {
		push @{$stock->gf_DE}, $data;
		warn "[parsing alignment: $data]\n";

	    } elsif ($field eq 'APH') {
		push @{$stock->gf_NH}, $data;

	    } elsif ($field eq 'NAM') {
		$name = $data;
		push @{$stock->seqname}, $name;
		$stock->seqdata->{$name} = "." x ($amino ? $aln : $dln);

	    } elsif ($field eq 'DSQ') {
		warn "$name: length is ", length($data), ", expected $dln"
		    if length($data) != $dln;
		$stock->seqdata->{$name} = $data
		    unless $amino;

	    } elsif ($field eq 'ASQ') {
		warn "$name: length is ", length($data), ", expected $aln"
		    if length($data) != $aln;
		if ($amino) {
		    $stock->seqdata->{$name} = $data;
		} else {
		    $stock->gr_AA->{$name} = amino1to3 ($data);
		}

	    } else {
#		warn "Unknown field '$field'\n";
		push @{$stock->gf_CC}, $line
		    unless $slim;

		if ($field eq 'ALN') { $aln = $data }
		if ($field eq 'DLN') { $dln = $data }
	    }

	} else {
	    warn "Ignoring line $line\n";
	}
    }

    my $lcols = $stock->lcols;
    my $maxcols = 80;
    while (($maxcols - 1 - $lcols) % 3 != 0) { --$maxcols }
    print $stock->to_string ($maxcols);
}


# convert single-letter into three-letter amino acid codes
sub amino1to3 {
    my ($seq) = @_;

    my %triplet = (A => Ala,
		   B => Asx,
		   C => Cys,
		   D => Asp,
		   E => Glu,
		   F => Phe,
		   G => Gly,
		   H => His,
		   I => Ile,
		   K => Lys,
		   L => Leu,
		   M => Met,
		   N => Asn,
		   P => Pro,
		   Q => Gln,
		   R => Arg,
		   S => Ser,
		   T => Thr,
		   V => Val,
		   W => Trp,
		   Y => Tyr
		   );

    my $ret = "";
    for (my $i = 0; $i < length($seq); ++$i) {
	my $c = uc substr ($seq, $i, 1);
	$ret .= exists ($triplet{$c}) ? $triplet{$c} : $c x 3;
    }

    return $ret;
}
